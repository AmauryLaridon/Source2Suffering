#%% --------------------------------------------------------------------------------#
# update for publication combined plot showing absolute cohort sizes and pie charts #
# ----------------------------------------------------------------------------------#

def plot_combined_population_update(
    df_GMT_strj,
    ds_pf_gs,
    da_gs_popdenom,
    gdf_country_borders,
    sims_per_step,
    flags,
    df_countries,
):

    # plot characteristics
    plot_var='unprec_99.99'
    x=12
    y=7
    markersize=10
    col_cbticlbl = 'gray'   # colorbar color of tick labels
    col_cbtic = 'gray'   # colorbar color of ticks
    col_cbedg = 'gray'   # colorbar color of edge
    cb_ticlen = 3.5   # colorbar length of ticks
    cb_ticwid = 0.4   # colorbar thickness of ticks
    cb_edgthic = 0   # colorbar thickness of edges between colors   
    by=2020
    cmap_whole = plt.cm.get_cmap('Reds')
    levels_cmap = np.arange(0,1.01,0.05)
    colors = [cmap_whole(i) for i in levels_cmap[:-1]]
    cmap_list_frac = mpl.colors.ListedColormap(colors,N=len(colors))
    ticks = np.arange(0,101,10)
    levels = np.arange(0,101,5)
    norm=mpl.colors.BoundaryNorm(levels,ncolors=len(levels)-1)
    l = 0 # letter indexing
    gmt_indices_152535 = [20,10,0]
    map_letters = {20:'c',10:'d',0:'e'}

    gmt_legend={
        GMT_indices_plot[0]:'1.5',
        GMT_indices_plot[1]:'2.5',
        GMT_indices_plot[2]:'3.5',
    }     

    f = plt.figure(figsize=(x,y))    
    gs0 = gridspec.GridSpec(10,3)

    # box plots
    # ax0 = f.add_subplot(gs0[0:4,0:2]) 
    ax0 = f.add_subplot(gs0[0:5,0:2]) 

    # pop totals
    # ax1 = f.add_subplot(gs0[4:,0:2],sharex=ax0)
    ax1 = f.add_subplot(gs0[5:,0:2],sharex=ax0)

    gs0.update(hspace=0)

    # maps
    gs01 = gridspec.GridSpecFromSubplotSpec(
        3,
        1, 
        subplot_spec=gs0[0:9,2:],
    )
    ax01 = f.add_subplot(gs01[0],projection=ccrs.Robinson())
    ax11 = f.add_subplot(gs01[1],projection=ccrs.Robinson())
    ax21 = f.add_subplot(gs01[2],projection=ccrs.Robinson()) 
    
    # this is copied from my pyramid plot maps. comparison meant to get nicer version of robinson seen in the other
    # seems that when robinson is called from subplot instantiation, the panels look different
    # f,ax = plt.subplots(
    #     ncols=1,
    #     nrows=1,
    #     subplot_kw={'projection':ccrs.Robinson()},
    #     transform=ccrs.PlateCarree()
    # )
    # ax.add_feature(feature.NaturalEarthFeature('physical', 'ocean', '50m', edgecolor='powderblue', facecolor='powderblue'))
    # gdf_p.to_crs(robinson).plot(
    #     ax=ax,
    #     column='grdi_q_by_p',
    #     color='darkgoldenrod',
    #     zorder=5,
    #     markersize=0.5,
    # )    
    # gdf_r.to_crs(robinson).plot(
    #     ax=ax,
    #     column='grdi_q_by_p',
    #     color='forestgreen',
    #     zorder=5,
    #     markersize=0.5,
    # )            
    # ax.set_xlim(gdf_robinson_bounds[0],gdf_robinson_bounds[2])
    # ax.set_ylim(gdf_robinson_bounds[1],gdf_robinson_bounds[3])      
    
    pos00 = ax21.get_position()
    cax00 = f.add_axes([
        pos00.x0-0.05,
        pos00.y0-0.1,
        pos00.width*1.95,
        pos00.height*0.2
    ])

    # get data
    df_list_gs = []
    extr='heatwavedarea'
    with open(data_dir+'{}/{}/isimip_metadata_{}_{}_{}.pkl'.format(flags['version'],extr,extr,flags['gmt'],flags['rm']), 'rb') as file:
        d_isimip_meta = pk.load(file)              
    with open(data_dir+'{}/{}/gridscale_aggregated_pop_frac_{}.pkl'.format(flags['version'],extr,extr), 'rb') as file:
        ds_pf_gs_plot = pk.load(file)
    da_p_gs_plot = ds_pf_gs_plot[plot_var].loc[{
        'GMT':GMT_indices_plot,
        'birth_year':sample_birth_years,
    }]
    sims_per_step = {}
    for step in GMT_labels:
        sims_per_step[step] = []
        for i in list(d_isimip_meta.keys()):
            if d_isimip_meta[i]['GMT_strj_valid'][step]:
                sims_per_step[step].append(i)  
                
    # --------------------------------------------------------------------
    # pf boxplot time series
    for step in GMT_indices_plot:
        da_pf_gs_plot_step = da_p_gs_plot.loc[{'run':sims_per_step[step],'GMT':step}].fillna(0).sum(dim='country') / da_gs_popdenom.sum(dim='country') * 100
        df_pf_gs_plot_step = da_pf_gs_plot_step.to_dataframe(name='pf').reset_index()
        df_pf_gs_plot_step['GMT_label'] = df_pf_gs_plot_step['GMT'].map(gmt_legend)       
        df_pf_gs_plot_step['hazard'] = extr
        df_list_gs.append(df_pf_gs_plot_step)
    df_pf_gs_plot = pd.concat(df_list_gs)

    
    colors = dict(zip(list(gmt_legend.values()),['steelblue','darkgoldenrod','darkred']))
    p = sns.boxplot(
        data=df_pf_gs_plot[df_pf_gs_plot['hazard']==extr],
        x='birth_year',
        y='pf',
        hue='GMT_label',
        palette=colors,
        showcaps=False,
        showfliers=False,
        boxprops={
            'linewidth':0,
            'alpha':0.5
        },        
        ax=ax0,
    )

    # Update #
    median_pf_per_birthyear = df_pf_gs_plot.groupby(["birth_year", "GMT_label"])["pf"].median()

    p.legend_.remove()                  
    ax0.spines['right'].set_visible(False)
    ax0.spines['top'].set_visible(False)      
    ax0.tick_params(colors='gray')
    ax0.set_ylim(0,100)
    ax0.spines['left'].set_color('gray')
    ax0.spines['bottom'].set_color('gray')      
    ax0.set_ylabel('$\mathregular{CF_{Heatwaves}}$ [%]',color='gray',fontsize=14)
    ax0.set_xlabel('Birth year',color='gray',fontsize=14)       
    ax0.set_title(
        letters[l],
        loc='left',
        fontweight='bold',
        fontsize=10
    )    
    l+=1 

    # bbox
    x0 = 0.025
    y0 = 0.7
    xlen = 0.2
    ylen = 0.3

    # space between entries
    legend_entrypad = 0.5

    # length per entry
    legend_entrylen = 0.75

    legend_font = 10
    legend_lw=3.5   

    legendcols = list(colors.values())
    handles = [
        Rectangle((0,0),1,1,color=legendcols[0]),\
        Rectangle((0,0),1,1,color=legendcols[1]),\
        Rectangle((0,0),1,1,color=legendcols[2])
    ]

    labels= [
        '1.5 °C GMT warming by 2100',
        '2.5 °C GMT warming by 2100',
        '3.5 °C GMT warming by 2100',    
    ]

    ax0.legend(
        handles, 
        labels, 
        bbox_to_anchor=(x0, y0, xlen, ylen), 
        loc = 'upper left',
        ncol=1,
        fontsize=legend_font, 
        labelcolor='gray',
        mode="expand", 
        borderaxespad=0.,\
        frameon=False, 
        columnspacing=0.05, 
        handlelength=legend_entrylen, 
        handletextpad=legend_entrypad
    )   

    # --------------------------------------------------------------------
    # populations

    ax1.spines['right'].set_visible(False)
    # ax1.spines['top'].set_visible(False)      
    ax1.tick_params(colors='gray')
    # ax1.set_ylim(0,100)
    ax1.spines['top'].set_color('gray')   
    ax1.spines['left'].set_color('gray')
    ax1.spines['bottom'].set_color('gray')    

    incomegroups = df_countries['incomegroup'].unique()
    income_countries = {}
    for category in incomegroups:
        income_countries[category] = list(df_countries.index[df_countries['incomegroup']==category])

    heights={}
    for category in incomegroups:
        heights[category] = da_gs_popdenom.loc[{
            'birth_year':np.arange(1960,2021,10),'country':income_countries[category]
        }].sum(dim='country').values / 10**6
    testdf = pd.DataFrame(heights)    
    testdf['birth_year'] = np.arange(1960,2021,10)
    testdf = testdf.set_index('birth_year')             
    testdf['total'] = testdf.sum(axis=1)
    p1 = testdf['total'].plot(
        kind='bar',
        # column='total',
        # stacked=True,
        color='gray',      
        ax=ax1,
        legend=False,
        rot=0.5
    ) 

    # Recover the medians for each GMT scenario
    median_1_5 = median_pf_per_birthyear.xs('1.5', level='GMT_label')
    median_2_5 = median_pf_per_birthyear.xs('2.5', level='GMT_label')
    median_3_5 = median_pf_per_birthyear.xs('3.5', level='GMT_label')

    # Creates the new bar plots for each scenario
    # Compute the total size of the global cohort that will be affected by scenario
    cohort_affected_1_5 = testdf['total'] * median_1_5.values / 100
    cohort_affected_2_5 = testdf['total'] * median_2_5.values / 100
    cohort_affected_3_5 = testdf['total'] * median_3_5.values / 100

    # Plots of the new bars per scenario

    p1_3 = cohort_affected_3_5.plot(
        kind='bar',
        # column='total',
        # stacked=True,
        color='darkred',      
        ax=ax1,
        legend=False,
        rot=0.5,
        stacked=True,
        bottom=0
    )

    p1_2 = cohort_affected_2_5.plot(
        kind='bar',
        # column='total',
        # stacked=True,
        color='darkgoldenrod',      
        ax=ax1,
        legend=False,
        rot=0.5,
        stacked=True
    )

    p1_1 = cohort_affected_1_5.plot(
        kind='bar',
        # column='total',
        # stacked=True,
        color='steelblue',      
        ax=ax1,
        legend=False,
        rot=0.5,
        stacked=True
    )

    #ax1.invert_yaxis()
    ax1.set_xlabel('Birth year',color='gray',fontsize=14)
    ax1.set_ylabel('Global cohort \n totals [in millions]',color='gray',fontsize=14)  
    ax1.set_title(
        letters[l],
        loc='left',
        fontweight='bold',
        fontsize=10
    )  
    # plot total cohort sizes as text inside bars
    for i,by in enumerate(testdf.index.values):
        ax1.text(
            x=i,
            y=np.round(testdf['total'].loc[by]) - 4.5,
            s=str(int(np.round(testdf['total'].loc[by]))),
            horizontalalignment='center',
            verticalalignment='center',
            # transform=ax1.transData,
            fontsize=8,
            color='white',
        )
    # plot total cohort sizes affected by the strongest impact scenario as text inside bars
    for i,by in enumerate(testdf.index.values):

        if i==0:
            ax1.text(
                x=i,
                y=cohort_affected_1_5.iloc[i] - 3.4,
                s=str(int(np.round(cohort_affected_3_5.iloc[i]))),
                horizontalalignment='center',
                verticalalignment='center',
                # transform=ax1.transData,
                fontsize=8,
                color='white',
            )
        elif i==1:
            ax1.text(
                x=i,
                y=cohort_affected_1_5.iloc[i] - 3.4,
                s=str(int(np.round(cohort_affected_1_5.iloc[i]))),
                horizontalalignment='center',
                verticalalignment='center',
                # transform=ax1.transData,
                fontsize=8,
                color='white',
            )
        elif i>1:
            ax1.text(
                x=i,
                y=cohort_affected_3_5.iloc[i] - 3.4,
                s=str(int(np.round(cohort_affected_3_5.iloc[i]))),
                horizontalalignment='center',
                verticalalignment='center',
                # transform=ax1.transData,
                fontsize=8,
                color='white',
            )

    # plot total cohort sizes affected by the strongest impact scenario as text inside bars
    for i,by in enumerate(testdf.index.values):

        if i<=1:
            pass
        else:
            ax1.text(
                x=i,
                y=cohort_affected_2_5.iloc[i] - 3.4,
                s=str(int(np.round(cohort_affected_2_5.iloc[i]))),
                horizontalalignment='center',
                verticalalignment='center',
                # transform=ax1.transData,
                fontsize=8,
                color='white',
            )
    
    # plot total cohort sizes affected by the strongest impact scenario as text inside bars
    for i,by in enumerate(testdf.index.values):

        if i<=1:
            pass
        else:
            ax1.text(
                x=i,
                y=cohort_affected_1_5.iloc[i] - 3.2,
                s=str(int(np.round(cohort_affected_1_5.iloc[i]))),
                horizontalalignment='center',
                verticalalignment='center',
                # transform=ax1.transData,
                fontsize=8,
                color='white',
            )
    
    
    l+=1
    # --------------------------------------------------------------------
    # maps of pop frac emergence for countries at 1, 2 and 3 deg pathways

    da_p_gs_plot = ds_pf_gs[plot_var].loc[{
        'GMT':gmt_indices_152535,
        'birth_year':by,
    }]
    df_list_gs = []
    for step in gmt_indices_152535:
        da_p_gs_plot_step = da_p_gs_plot.loc[{'run':sims_per_step[step],'GMT':step}].mean(dim='run')
        da_p_gs_plot_step = da_p_gs_plot_step / da_gs_popdenom.loc[{'birth_year':by}] * 100
        df_p_gs_plot_step = da_p_gs_plot_step.to_dataframe(name='pf').reset_index()
        df_p_gs_plot_step = df_p_gs_plot_step.assign(GMT_label = lambda x: np.round(df_GMT_strj.loc[2100,x['GMT']],1).values.astype('str'))
        df_list_gs.append(df_p_gs_plot_step)
    df_p_gs_plot = pd.concat(df_list_gs)
    df_p_gs_plot['pf'] = df_p_gs_plot['pf'].fillna(0)  
    gdf = cp(gdf_country_borders.reset_index())
    gdf_p = cp(gdf_country_borders.reset_index())
    robinson = ccrs.Robinson().proj4_init

    for ax,step in zip((ax01,ax11,ax21),gmt_indices_152535):
        gdf_p['pf']=df_p_gs_plot['pf'][df_p_gs_plot['GMT']==step].values
        ax.add_feature(feature.NaturalEarthFeature('physical', 'ocean', '50m', edgecolor='powderblue', facecolor='powderblue'))
        gdf_p.to_crs(robinson).plot(
            ax=ax,
            column='pf',
            cmap=cmap_list_frac,
            norm=norm,
            cax=cax00,
            zorder=2,
            rasterized=True,
        )

        gdf.to_crs(robinson).plot(
            ax=ax,
            color='none', 
            edgecolor='black',
            linewidth=0.25,
            zorder=3,
        )
        
        gdf_robinson_bounds = gdf_p.to_crs(robinson).total_bounds # (minx,miny,maxx,maxy)
        
        # ax.set_global() # checking to see if this makes projections fully showing antarctica
        
        ax.set_title(
            '{} °C'.format(gmt_legend[step]),
            loc='center',
            fontweight='bold',
            fontsize=12,
            color='gray',       
        )
        
        ax.set_title(
            map_letters[step],
            loc='left',
            fontweight='bold',
            fontsize=10
        )    

        # triangles showing connections between GMT step box plot to map panels ---------------
        from matplotlib.patches import Circle, Wedge, Polygon
        from matplotlib.collections import PatchCollection
        if step == gmt_indices_152535[0]:
            x_h=1 
        elif step == gmt_indices_152535[1]:
            x_h=0.95                      
        elif step == gmt_indices_152535[-1]:
            x_h=0.9
        y_h= df_pf_gs_plot[(df_pf_gs_plot['birth_year']==by)&(df_pf_gs_plot['GMT']==step)]['pf'].median() / 100
        x_m=0
        y_m_i=0
        y_m_f=1.0
        con_low = ConnectionPatch(
            xyA=(x_h,y_h),
            xyB=(x_m,y_m_i),
            coordsA=ax0.transAxes,
            coordsB=ax.transAxes,
            color='lightgray',
            alpha=0.5,
            zorder=5,
        )
        con_hi = ConnectionPatch(
            xyA=(x_h,y_h),
            xyB=(x_m,y_m_f),
            coordsA=ax0.transAxes,
            coordsB=ax.transAxes,
            color='lightgray',
            alpha=0.5,
            zorder=5,
        )   
        con_vert = ConnectionPatch(
            xyA=(x_m,y_m_i),
            xyB=(x_m,y_m_f),
            coordsA=ax.transAxes,
            coordsB=ax.transAxes,
            color='lightgray',
            alpha=0.5,
            zorder=5,
        )
        
        line_low = con_low.get_path().vertices
        line_hi = con_hi.get_path().vertices
        
        ax0.add_artist(con_low)  
        ax0.add_artist(con_hi) 
        
        tri_coords = np.stack(
            (con_low.get_path().vertices[0],con_low.get_path().vertices[-1],con_hi.get_path().vertices[-1]),
            axis=0,
        )
        triangle = plt.Polygon(tri_coords,ec='lightgray',fc='lightgray',alpha=0.5,zorder=10,clip_on=False)    
        ax0.add_artist(triangle)           
        
        # triangles showing connections between GMT step box plot to map panels ---------------        
        
    cb = mpl.colorbar.ColorbarBase(
        ax=cax00, 
        cmap=cmap_list_frac,
        norm=norm,
        orientation='horizontal',
        spacing='uniform',
        ticks=ticks,
        drawedges=False,
    )

    cb.set_label(
        '$\mathregular{CF_{Heatwaves}}$ for 2020 birth cohort [%]',
        fontsize=14,
        color='gray'
    )
    cb.ax.xaxis.set_label_position('top')
    cb.ax.tick_params(
        labelcolor=col_cbticlbl,
        color=col_cbtic,
        length=cb_ticlen,
        width=cb_ticwid,
        direction='out'
    )   
    cb.outline.set_edgecolor(col_cbedg)
    cb.outline.set_linewidth(cb_edgthic)   
    cax00.xaxis.set_label_position('top')   

    f.savefig(scripts_dir+'/figures/grant_2025/ms_figures/f2_combined_plot_popsizes_update.png',dpi=1000)
    f.savefig(scripts_dir+'/figures/grant_2025/ms_figures/f2_combined_plot_popsizes_update.pdf', format="pdf", bbox_inches="tight")
    return gdf_robinson_bounds