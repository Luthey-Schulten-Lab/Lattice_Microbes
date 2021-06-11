# wget https://developer.nvidia.com/compute/cuda/9.2/Prod/local_installers/cuda_9.2.88_396.26_linux -O /tmp/cudatoolkit.sh &&
# wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O /tmp/miniconda.sh &&

userCfg = <<-SHELL
echo "export PATH=/usr/local/cuda/bin:/home/vagrant/miniconda/bin:/usr/local/bin:\$PATH" >> /home/vagrant/.bashrc &&
echo "source /home/vagrant/miniconda/bin/activate root" >> /home/vagrant/.bashrc &&
bash /tmp/miniconda.sh -b -p /home/vagrant/miniconda &&
rm /tmp/miniconda.sh &&
/home/vagrant/miniconda/bin/conda install -y protobuf hdf5 numpy cython cmake ninja swig &&
/home/vagrant/miniconda/bin/conda install -y -c anaconda git gcc_linux-64 gxx_linux-64 binutils_linux-64
SHELL

sysCfg =  <<-SHELL
bash /tmp/cudatoolkit.sh --silent --toolkit --override &&
rm /tmp/cudatoolkit.sh
SHELL

cudaFile = "vagrant-dist/cuda_9.2.88_396.26_linux"
condaFile = "vagrant-dist/Miniconda3-latest-Linux-x86_64.sh"

Vagrant.configure("2") do |config|
    config.vm.define "centos5" do |centos5|
        centos5.vm.box = "bento/centos-5.11"
        centos5.vm.provision "file", source: cudaFile, destination: "/tmp/cudatoolkit.sh"
        centos5.vm.provision "file", source: condaFile, destination: "/tmp/miniconda.sh"
        centos5.vm.provision "shell", inline: sysCfg
        centos5.vm.provision "shell", privileged: false, inline: userCfg
        centos5.vm.synced_folder "src", "/vagrant"
    end
    config.vm.define "centos6" do |centos6|
        centos6.vm.box = "centos/6"
        centos6.vm.provision "file", source: cudaFile, destination: "/tmp/cudatoolkit.sh"
        centos6.vm.provision "file", source: condaFile, destination: "/tmp/miniconda.sh"
        centos6.vm.provision "shell", inline: sysCfg
        centos6.vm.provision "shell", privileged: false, inline: userCfg
        centos6.vm.synced_folder "src", "/vagrant", type: "rsync", rsync__exclude: ".git/"
    end
    config.vm.define "centos7" do |centos7|
        centos7.vm.box = "centos/7"
        centos7.vm.provision "file", source: cudaFile, destination: "/tmp/cudatoolkit.sh"
        centos7.vm.provision "file", source: condaFile, destination: "/tmp/miniconda.sh"
        centos7.vm.provision "shell", inline: sysCfg
        centos7.vm.provision "shell", privileged: false, inline: userCfg
        centos7.vm.synced_folder "src", "/vagrant", type: "rsync", rsync__exclude: ".git/"
    end
    config.vm.define "ubuntu-trusty" do |trusty|
        trusty.vm.box = "ubuntu/trusty64"
        trusty.vm.provision "file", source: cudaFile, destination: "/tmp/cudatoolkit.sh"
        trusty.vm.provision "file", source: condaFile, destination: "/tmp/miniconda.sh"
        trusty.vm.provision "shell", inline: sysCfg
        trusty.vm.provision "shell", privileged: false, inline: userCfg
        trusty.vm.synced_folder "src", "/vagrant", type: "rsync", rsync__exclude: ".git/"
    end
    config.vm.define "ubuntu-xenial" do |xenial|
        xenial.vm.box = "ubuntu/xenial64"
        xenial.vm.provision "file", source: cudaFile, destination: "/tmp/cudatoolkit.sh"
        xenial.vm.provision "file", source: condaFile, destination: "/tmp/miniconda.sh"
        xenial.vm.provision "shell", inline: sysCfg
        xenial.vm.provision "shell", privileged: false, inline: userCfg
        xenial.vm.synced_folder "src", "/vagrant", type: "rsync", rsync__exclude: ".git/"
    end
    config.vm.define "ubuntu-precise" do |precise|
        precise.vm.box = "ubuntu/precise64"
        precise.vm.provision "file", source: cudaFile, destination: "/tmp/cudatoolkit.sh"
        precise.vm.provision "file", source: condaFile, destination: "/tmp/miniconda.sh"
        precise.vm.provision "shell", inline: sysCfg
        precise.vm.provision "shell", privileged: false, inline: userCfg
        precise.vm.synced_folder "src", "/vagrant", type: "rsync", rsync__exclude: ".git/"
    end
end
