function exportgraphics19(filename)
    % Export current figure to PDF (obsolete solution for MATLAB R2019b)
    fighandle=gcf;
    fighandle.PaperOrientation='landscape';
    %fighandle.PaperType='uslegal';
    fighandle.PaperSize=[18 8];
    print(filename,'-dpdf','-fillpage');
    system(['pdfcrop ' filename '.pdf ' filename '.pdf'])