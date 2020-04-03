require(multipanelfigure)

fig <- multi_panel_figure(columns = 4, rows = 2, width = 297, height = 210)

filename <- '../multicrisprout/A_crispr.pdf'
file.copy('Z:/PaperInPrep/multicrispr/figures/fig1/A_crispr.pdf', filename)
image <- magick::image_read_pdf(filename)
filename %<>% stringi::stri_replace_first_fixed('pdf', 'png')
magick::image_write(image, filename)
fig <- fill_panel(fig, filename, row = 1, column = 1:2, scaling = 'shrink')

filename <- '../multicrisprout/B_prime_editing.pdf'
file.copy('Z:/PaperInPrep/multicrispr/figures/fig1/B_prime_editing.pdf', filename)
image <- magick::image_read_pdf(filename)
filename %<>% stringi::stri_replace_first_fixed('pdf', 'png')
magick::image_write(image, filename)
fig <- fill_panel(fig, filename, row = 1, column = 3:4, scaling = 'shrink')
fig

filename <- '../multicrisprout/C_grna_tools.pdf'
file.copy('Z:/PaperInPrep/multicrispr/figures/fig1/C_grna_tools.pdf', filename)
image <- magick::image_read_pdf(filename)
filename %<>% stringi::stri_replace_first_fixed('pdf', 'png')
magick::image_write(image, filename)
fig <- fill_panel(fig, filename, row = 2, column = 1, scaling = 'shrink')
fig

filename <- '../multicrisprout/D_genome_arithmetics_wide.pdf'
file.copy('Z:/PaperInPrep/multicrispr/figures/fig1/D_genome_arithmetics_wide.pdf', filename)
image <- magick::image_read_pdf(filename)
filename %<>% stringi::stri_replace_first_fixed('pdf', 'png')
magick::image_write(image, filename)
fig <- fill_panel(fig, filename, row = 2, column = 2:4, scaling = 'shrink')
fig

pdf('Z:/PaperInPrep/multicrispr/figures/fig1/fig1.pdf', 
    width  = figure_width( fig, 'inch'), 
    height = figure_height(fig, 'inch'))
print(fig)
dev.off()
