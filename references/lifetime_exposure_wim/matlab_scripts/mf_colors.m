

% --------------------------------------------------------------------
% function to define some nice colors
% --------------------------------------------------------------------


function [colors] = mf_colors


% nice hex combinations are taken from:
% https://designschool.canva.com/blog/website-color-schemes/
% https://designschool.canva.com/blog/100-color-combinations/


% line colors
colors    = [0.20 0.73 0.03     ; ...   % 1.  grass 
             0.82 0.41 0.12     ; ...   % 2.  brown ('CHOCOLATE' in mf_rgb)
             0.80 0.73 0.03     ; ...   % 3.  sand
             1.00 0.00 0.00     ; ...   % 4.  red
             0.60 0.00 0.40     ; ...   % 5.  purple
             1.00 0.40 1.00     ; ...   % 6.  pink
             1.00 0.80 0.00     ; ...   % 7.  yellow
             0.50 0.50 0.50     ; ...   % 8.  grey
             0.00 0.00 1.00     ; ...   % 9.  blue
             0.30 0.00 0.00     ; ...   % 10. brown
             1.00 0.40 0.10     ; ...   % 11. orange
             0.75 0.75 1.00     ; ...   % 12. light blue
             1.00 0.75 0.75     ; ...   % 13. light red
             0.75 1.00 0.75     ; ...   % 14. light green
             0.90 0.83 0.13     ; ...   % 15. light sand
             0.89 0.10 0.11     ; ...   % 16. CB red
             0.21 0.49 0.72     ; ...   % 17. CB blue
             0.30 0.68 0.29     ; ...   % 18. CB green
             0.99 0.73 0.52     ; ...   % 19. CB lightred
             0.65 0.74 0.86     ; ...   % 20. CB lightblue
             0.70 0.87 0.54     ; ...   % 21. CB lightgreen
             0.60 0.24 0.24     ; ...   % 22. nice red   (hsv2rgb([0   0.6 0.6]))
             0.31 0.60 0.24     ; ...   % 23. nice green (hsv2rgb([0.3 0.6 0.6]))
             0.24 0.38 0.60     ; ...   % 24. nice blue  (hsv2rgb([0.6 0.6 0.6]))
             
             hex2rgb('#2E4600') ; ...   % 25. Forest Green (combo 6)
             hex2rgb('#486B00') ; ...   % 26. Grass (combo 6)
             hex2rgb('#A2C523') ; ...   % 27. lime (combo 6)
             hex2rgb('#7D4427') ; ...   % 28. earth (combo 6)
             
             hex2rgb('#000000') ; ...   % 29. Black (combo 7)
             hex2rgb('#062F4F') ; ...   % 30. Ink (combo 7)
             hex2rgb('#813772') ; ...   % 31. Posy (combo 7)
             hex2rgb('#BB2601') ; ...   % 32. Embers (combo 7)
             
             hex2rgb('#4C3F54') ; ...   % 33. Fig (combo 36)
             hex2rgb('#D13525') ; ...   % 34. Apple red (combo 36)
             hex2rgb('#F2C057') ; ...   % 35. Swiss cheese (combo 36)
             hex2rgb('#486824') ; ...   % 36. Basil (combo 36)
             
             hex2rgb('#444C5C') ; ...   % 37. Faded Navy (combo 0)
             hex2rgb('#CE5A57') ; ...   % 38. Punch (combo 0)
             hex2rgb('#78A5A3') ; ...   % 39. Ocean Breeze (combo 0)
             hex2rgb('#E1B16A') ; ...   % 40. Warm (combo 0)
             
             hex2rgb('#003b46') ; ...   % 41. Deep Aqua (combo 5 - Cool Blues)
             hex2rgb('#07575b') ; ...   % 42. Ocean (combo 5 - Cool Blues)
             hex2rgb('#66a5ad') ; ...   % 43. Wave (combo 5 - Cool Blues)
             hex2rgb('#c4dfe6') ; ...   % 44. Seafoam (combo 5 - Cool Blues)

             hex2rgb('#46211a') ; ...   % 45. Crevice (combo 3 - Dark & Earthy)
             hex2rgb('#693d3d') ; ...   % 46. Cloud Shadow (combo 3 - Dark & Earthy)
             hex2rgb('#ba5536') ; ...   % 47. Desert (combo 3 - Dark & Earthy)
             hex2rgb('#a43820') ; ...   % 48. Red Clay (combo 3 - Dark & Earthy)

             hex2rgb('#8c0004') ; ...   % 49. Currant (combo 71 - Mediterranean afternoon)
             hex2rgb('#c8000a') ; ...   % 50. Scarlet (combo 71 - Mediterranean afternoon)
             hex2rgb('#e8a735') ; ...   % 51. Marigold (combo 71 - Mediterranean afternoon)
             hex2rgb('#e2c499') ; ...   % 52. Cobblestone (combo 71 - Mediterranean afternoon)

             0.02 0.62 0.47     ; ...   % 53. nice green (fig. 3a in Zekollari et al., 2020)
             0.46 0.44 0.70     ; ...   % 54. nice purple (fig. 3a in Zekollari et al., 202)             
                                        ]; 


end

