switch CHCase
    case 1
        [u,v]  =D2InitialUVza(U0,u,v);  %Zalesak's Disk时开启,该算例通用性不强，参考改写E:\Wangxs\Multiphase Flow\phase move CH\phase move CH2 zalesak
        [u] = D2set_BCNeu(u);[v] = D2set_BCNeu(v);
        uf=D2fXmean(u);vf=D2fYmean(v);
        uL=u;vL=v;
        ufL=uf ;  vfL=vf;  %插值面心速度
end
        phi_initial=phi;%用于定量对比