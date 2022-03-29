PRO SHARP_MOVIE,segment,files_dir,save_dir

;+
;
; Name SHARP_MOVIE
;      
; Purpose:
;          Takes a directory and makes a movie of a sharps segment file
;          with scc_mkmovie so the final file is .mvi and has the coordinates
;          and can be load with scc_playmoviem
;      
; Calling sequence:
;
; Input:
;          segment:  string type of segment to use for the movie 
;                    It can be one of several tipes
;                          Br, Bp, Bt, bitmap    
;                      
;
;          files_dir: directory where files are at
;          save_dir:  directory where movie is to be saved
;     
; Output:
;        .mvi movie file to be saved in the save_dir 
;
; Keywords:
;
;
; Author: Merenda C. Luciano A. 
;
;-

sharp_files = file_search(files_dir + '/*.' + segment + '.fits')

harpnum = sharp_files[16:20]

read_sdo,sharp_files[0],hdr,img

SCC_MKMOVIE,sharp_files,min(img),max(img),'hmi',/TIMES,SDIR=save_dir,SAVE=''

