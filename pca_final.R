#m=read.table("/Users/ds/columbia/StatUniverse/final_project/R-spect/SSGSS_info-R.txt",head=T)
 m=read.table("/home/rubab/geco/pca/infra_red_data/SSGSS_info-R.txt",head=T)

spect=matrix(0,331,97)
error=matrix(0,331,97)
lambda=c(0:330)*.1+5

for (i in c(1:97))
{
#filename<-paste("/Users/ds/columbia/StatUniverse/final_project/R-spect/SSGSS_",format(i),"-R.txt",sep="");
     filename<-paste("/home/rubab/geco/pca/infra_red_data/SSGSS_",format(i),"-R.txt",sep="");
     d=read.table(filename,head=T);

     sp_int=approx(d$wave/(1+m$z[i]),d$flux,lambda)
     spect[1:331,i]=sp_int$y

     err_int=approx(d$wave/(1+m$z[i]),d$error,lambda)
     error[1:331,i]=err_int$y
}

median_spect=matrix(0,331,1)
for (i in c(1:331))
{
median_spect[i,1]=median(spect[i,],na.rm=T);
}

median_matrix=matrix(0,331,97)
for (i in c(1:97))
{
median_matrix[1:331,i]=median_spect[1:331]
}

adj_spect=matrix(0,300,97)
adj_spect = (spect[6:305,] - median_matrix[6:305,])

nlines = 21
line_wavelength = c(1:nlines)*0
line_wavelength[1]  = 6.2    # PAH        PAH
line_wavelength[2]  = 7.0    # [ArII]     FS 
line_wavelength[3]  = 7.7    # PAH        PAH
line_wavelength[4]  = 8.6    # PAH        PAH
line_wavelength[5]  = 9.0    # [ArIII]    FS
line_wavelength[6]  = 9.7    # H2         Mol.
line_wavelength[7]  = 10.5   # [SIV]      FS
line_wavelength[8]  = 11.3   # PAH        PAH
line_wavelength[9]  = 12.0   # PAH        PAH
line_wavelength[10] = 12.3   # H2         Mol.
line_wavelength[11] = 12.7   # PAH        PAH
line_wavelength[12] = 12.8   # [NeII]     FS
line_wavelength[13] = 15.6   # [NeIII]    FS
line_wavelength[14] = 16.5   # PAH        PAH
line_wavelength[15] = 17.0   # H2         Mol.
line_wavelength[16] = 17.9   # [FeII]     FS
line_wavelength[17] = 18.3   # [NeIII]    FS
line_wavelength[18] = 18.8   # [SIII]     FS
line_wavelength[19] = 26.0   # [OIV/FeII] FS
line_wavelength[20] = 33.5   # [SIII]     FS
line_wavelength[21] = 34.8   # [SiII]     FS


# *********** Method-A ************

k = 1
cov_matrix = matrix(0,300,300)
for (i in c(1:300))
{
	for (j in c(1:300))
	{
		cov_matrix[i,j] = cov(adj_spect[i,1:97],adj_spect[j,1:97],use="complete.obs")
	}
}

eigen_stuff = eigen(cov_matrix)
eigen_values = eigen_stuff$values

feature_vector_matrix = eigen_stuff$vectors
for (i in c(1:3))
{
	plot_name = sprintf("Method A: Eigen Spectra - %d",i)
	fig_name = sprintf("Figure - %d",k)
	file_name= sprintf("Figure_%d.eps",k)

	postscript(file=file_name,title=fig_name)
	plot(lambda[6:305],(feature_vector_matrix[1:300,i]*(-1)+median_spect[6:305]),ty="l",col="purple")

	for (i in c(1:nlines))
	{
		lines(c(line_wavelength[i],line_wavelength[i]),c(-100,100),lty=1,col="red")
	}

	title(main=fig_name,sub=plot_name)
	dev.off()
	k = k+1
}

transposed_feature_vector_matrix = t(feature_vector_matrix)
new_data = (transposed_feature_vector_matrix[1:3,]%*%adj_spect)
recovered_orig_data = (feature_vector_matrix[,1:3]%*%new_data)

for (i in c(1,10,25,50,97))
{
	plot_name = sprintf("Method A: Original (blue) and Recovered (green) Spectra of Galaxy - %d",i)
	fig_name = sprintf("Figure - %d",k)
	file_name= sprintf("Figure_%d.eps",k)

	postscript(file=file_name,title=fig_name)
	plot(lambda[6:305],(adj_spect[1:300,i]+median_spect[6:305]),ty="l",col="blue")
	lines(lambda[6:305],(recovered_orig_data[1:300,i]+median_spect[6:305]),ty="l",col="green")
	
	for (i in c(1:nlines))
	{
		lines(c(line_wavelength[i],line_wavelength[i]),c(-100,100),lty=1,col="red")
	}

	title(main=fig_name,sub=plot_name)
	dev.off()
	k = k+1
}


# *********** Method-B ************

cov_matrix = matrix(0,97,97)
for (i in c(1:97))
{
	for (j in c(1:97))
	{
		cov_matrix[i,j] = cov(adj_spect[1:300,i],adj_spect[1:300,j],use="complete.obs")
	}
}

eigen_stuff = eigen(cov_matrix)
eigen_values = eigen_stuff$values
feature_vector_matrix = eigen_stuff$vectors
transposed_feature_vector_matrix = t(feature_vector_matrix)
transposed_adj_spect = t(adj_spect)

new_data = (transposed_feature_vector_matrix[1:3,]%*%transposed_adj_spect)
transposed_new_data = t(new_data)
for (i in c(1:3))
{
	plot_name = sprintf("Method B: Eigen Spectra - %d",i)
	fig_name = sprintf("Figure - %d",k)
	file_name= sprintf("Figure_%d.eps",k)

	postscript(file=file_name,title=fig_name)
	plot(lambda[6:305],(transposed_new_data[1:300,i]+median_spect[6:305]),ty="l",col="purple")

	for (i in c(1:nlines))
	{
		lines(c(line_wavelength[i],line_wavelength[i]),c(-100,100),lty=1,col="red")
	}

	title(main=fig_name,sub=plot_name)
	dev.off()
	k = k+1
}

recovered_orig_data = (feature_vector_matrix[,1:3]%*%new_data)
transposed_recovered_orig_data = t(recovered_orig_data)

for (i in c(1,10,25,50,97))
{
	plot_name = sprintf("Method B: Original (blue) and Recovered (green) Spectra of Galaxy - %d",i)
	fig_name = sprintf("Figure - %d",k)
	file_name= sprintf("Figure_%d.eps",k)

	postscript(file=file_name,title=fig_name)
	plot(lambda[6:305],(adj_spect[1:300,i]+median_spect[6:305]),ty="l",col="blue")
	lines(lambda[6:305],(transposed_recovered_orig_data[1:300,i]+median_spect[6:305]),ty="l",col="green")
	
	for (i in c(1:nlines))
	{
		lines(c(line_wavelength[i],line_wavelength[i]),c(-100,100),lty=1,col="red")
	}

	title(main=fig_name,sub=plot_name)
	dev.off()
	k = k+1
}

