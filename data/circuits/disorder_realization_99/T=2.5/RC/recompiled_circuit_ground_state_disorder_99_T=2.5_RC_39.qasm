OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.34243256) q[0];
sx q[0];
rz(2.7437796) q[0];
sx q[0];
rz(7.5510511) q[0];
rz(-0.59029382) q[1];
sx q[1];
rz(2.3550912) q[1];
sx q[1];
rz(13.785706) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40400051) q[0];
sx q[0];
rz(-2.1426179) q[0];
sx q[0];
rz(0.19947995) q[0];
rz(-pi) q[1];
rz(-2.8706495) q[2];
sx q[2];
rz(-0.56519714) q[2];
sx q[2];
rz(0.2172367) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.089069033) q[1];
sx q[1];
rz(-0.34862754) q[1];
sx q[1];
rz(-2.2050956) q[1];
x q[2];
rz(-2.9079535) q[3];
sx q[3];
rz(-2.2147182) q[3];
sx q[3];
rz(-1.131191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.56403247) q[2];
sx q[2];
rz(-1.6200248) q[2];
sx q[2];
rz(1.3684028) q[2];
rz(-2.2300301) q[3];
sx q[3];
rz(-1.7715958) q[3];
sx q[3];
rz(-0.024356775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0381222) q[0];
sx q[0];
rz(-0.39762527) q[0];
sx q[0];
rz(0.11904112) q[0];
rz(-2.0114404) q[1];
sx q[1];
rz(-0.96769133) q[1];
sx q[1];
rz(-1.4281323) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85066909) q[0];
sx q[0];
rz(-0.68456283) q[0];
sx q[0];
rz(0.67770358) q[0];
rz(-pi) q[1];
x q[1];
rz(1.150564) q[2];
sx q[2];
rz(-1.9398076) q[2];
sx q[2];
rz(2.8431161) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.88682879) q[1];
sx q[1];
rz(-2.6568251) q[1];
sx q[1];
rz(-1.280275) q[1];
rz(-2.302243) q[3];
sx q[3];
rz(-0.4412868) q[3];
sx q[3];
rz(1.6755144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6156562) q[2];
sx q[2];
rz(-1.6776513) q[2];
sx q[2];
rz(0.76457912) q[2];
rz(0.48314759) q[3];
sx q[3];
rz(-1.8254447) q[3];
sx q[3];
rz(-2.29276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1044384) q[0];
sx q[0];
rz(-0.26532441) q[0];
sx q[0];
rz(1.6695439) q[0];
rz(-0.56023359) q[1];
sx q[1];
rz(-1.3872223) q[1];
sx q[1];
rz(-1.4405506) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9870696) q[0];
sx q[0];
rz(-1.8049585) q[0];
sx q[0];
rz(-2.9375141) q[0];
rz(-1.57651) q[2];
sx q[2];
rz(-2.2457321) q[2];
sx q[2];
rz(2.142148) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6145103) q[1];
sx q[1];
rz(-1.0497112) q[1];
sx q[1];
rz(0.39343843) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7974924) q[3];
sx q[3];
rz(-1.3406702) q[3];
sx q[3];
rz(3.074444) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.6560087) q[2];
sx q[2];
rz(-1.6501082) q[2];
sx q[2];
rz(0.34240016) q[2];
rz(-1.418669) q[3];
sx q[3];
rz(-0.72903052) q[3];
sx q[3];
rz(0.5595783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2320084) q[0];
sx q[0];
rz(-0.37264687) q[0];
sx q[0];
rz(1.5268071) q[0];
rz(-2.3341663) q[1];
sx q[1];
rz(-1.5734943) q[1];
sx q[1];
rz(1.9400914) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0394131) q[0];
sx q[0];
rz(-1.2713968) q[0];
sx q[0];
rz(0.79203006) q[0];
rz(-pi) q[1];
rz(2.3551201) q[2];
sx q[2];
rz(-2.3209089) q[2];
sx q[2];
rz(-0.72696668) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.11749841) q[1];
sx q[1];
rz(-1.9650241) q[1];
sx q[1];
rz(0.82347639) q[1];
rz(-0.19530794) q[3];
sx q[3];
rz(-0.26104078) q[3];
sx q[3];
rz(-1.9633499) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.80369192) q[2];
sx q[2];
rz(-0.98946977) q[2];
sx q[2];
rz(-1.9039512) q[2];
rz(-1.8303653) q[3];
sx q[3];
rz(-1.6236191) q[3];
sx q[3];
rz(-1.9633912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0230947) q[0];
sx q[0];
rz(-1.3730405) q[0];
sx q[0];
rz(-0.9088687) q[0];
rz(0.88242775) q[1];
sx q[1];
rz(-0.71417037) q[1];
sx q[1];
rz(0.73208255) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86504146) q[0];
sx q[0];
rz(-2.2892729) q[0];
sx q[0];
rz(-0.9263692) q[0];
x q[1];
rz(1.2901487) q[2];
sx q[2];
rz(-0.82260231) q[2];
sx q[2];
rz(1.6463239) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.27837023) q[1];
sx q[1];
rz(-1.178982) q[1];
sx q[1];
rz(-2.6696221) q[1];
x q[2];
rz(2.3084435) q[3];
sx q[3];
rz(-0.97104302) q[3];
sx q[3];
rz(2.7540327) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0714134) q[2];
sx q[2];
rz(-0.81255239) q[2];
sx q[2];
rz(-1.8522235) q[2];
rz(-1.6728801) q[3];
sx q[3];
rz(-1.2082992) q[3];
sx q[3];
rz(2.7220791) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5759204) q[0];
sx q[0];
rz(-1.9603632) q[0];
sx q[0];
rz(1.1015724) q[0];
rz(-2.9167602) q[1];
sx q[1];
rz(-1.79554) q[1];
sx q[1];
rz(0.98519957) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.12041872) q[0];
sx q[0];
rz(-0.93909953) q[0];
sx q[0];
rz(2.9178502) q[0];
x q[1];
rz(1.9198138) q[2];
sx q[2];
rz(-0.98528457) q[2];
sx q[2];
rz(-1.7988009) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.8781153) q[1];
sx q[1];
rz(-2.2190071) q[1];
sx q[1];
rz(-2.1489759) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.043126338) q[3];
sx q[3];
rz(-2.6131328) q[3];
sx q[3];
rz(0.915574) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.78099403) q[2];
sx q[2];
rz(-0.66086078) q[2];
sx q[2];
rz(-1.9471656) q[2];
rz(-3.0418975) q[3];
sx q[3];
rz(-1.7710779) q[3];
sx q[3];
rz(-2.1626507) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7128971) q[0];
sx q[0];
rz(-2.6850057) q[0];
sx q[0];
rz(-1.1140484) q[0];
rz(1.3975337) q[1];
sx q[1];
rz(-1.3797398) q[1];
sx q[1];
rz(0.31375113) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6614051) q[0];
sx q[0];
rz(-1.620694) q[0];
sx q[0];
rz(1.2779092) q[0];
rz(-pi) q[1];
rz(-0.6724486) q[2];
sx q[2];
rz(-0.20951665) q[2];
sx q[2];
rz(-2.5447951) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3849029) q[1];
sx q[1];
rz(-0.80763615) q[1];
sx q[1];
rz(2.7212423) q[1];
x q[2];
rz(-1.3082771) q[3];
sx q[3];
rz(-2.2394925) q[3];
sx q[3];
rz(1.4895542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.067001192) q[2];
sx q[2];
rz(-2.7374697) q[2];
sx q[2];
rz(2.5355549) q[2];
rz(1.0780942) q[3];
sx q[3];
rz(-2.2793016) q[3];
sx q[3];
rz(-1.3585453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.17230497) q[0];
sx q[0];
rz(-0.91738874) q[0];
sx q[0];
rz(2.0148) q[0];
rz(2.2276095) q[1];
sx q[1];
rz(-2.7028449) q[1];
sx q[1];
rz(-2.9235358) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.954531) q[0];
sx q[0];
rz(-2.0045223) q[0];
sx q[0];
rz(1.0591255) q[0];
rz(-pi) q[1];
rz(-0.86037068) q[2];
sx q[2];
rz(-2.2220774) q[2];
sx q[2];
rz(1.2597348) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.40849028) q[1];
sx q[1];
rz(-2.4839851) q[1];
sx q[1];
rz(1.3642743) q[1];
rz(3.0239437) q[3];
sx q[3];
rz(-1.7075286) q[3];
sx q[3];
rz(-2.4916388) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.5128936) q[2];
sx q[2];
rz(-0.67535526) q[2];
sx q[2];
rz(0.2891573) q[2];
rz(0.9946) q[3];
sx q[3];
rz(-1.1887487) q[3];
sx q[3];
rz(0.4579671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4299803) q[0];
sx q[0];
rz(-1.0884322) q[0];
sx q[0];
rz(-0.76643884) q[0];
rz(-1.537716) q[1];
sx q[1];
rz(-1.3130554) q[1];
sx q[1];
rz(2.2235353) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.799918) q[0];
sx q[0];
rz(-2.3324488) q[0];
sx q[0];
rz(2.2412712) q[0];
rz(-pi) q[1];
x q[1];
rz(0.16854281) q[2];
sx q[2];
rz(-1.4282116) q[2];
sx q[2];
rz(-3.0592205) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.2868216) q[1];
sx q[1];
rz(-0.76463759) q[1];
sx q[1];
rz(0.36061339) q[1];
rz(-pi) q[2];
rz(3.0357846) q[3];
sx q[3];
rz(-1.8195767) q[3];
sx q[3];
rz(0.63423587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.90091577) q[2];
sx q[2];
rz(-0.46611163) q[2];
sx q[2];
rz(-0.051699836) q[2];
rz(-0.85203552) q[3];
sx q[3];
rz(-1.7107191) q[3];
sx q[3];
rz(0.26028546) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.32605115) q[0];
sx q[0];
rz(-1.4651848) q[0];
sx q[0];
rz(-3.0503804) q[0];
rz(0.7971898) q[1];
sx q[1];
rz(-1.3858831) q[1];
sx q[1];
rz(-0.43112722) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2225616) q[0];
sx q[0];
rz(-3.0662144) q[0];
sx q[0];
rz(-2.4576709) q[0];
rz(-3.0431137) q[2];
sx q[2];
rz(-2.2171671) q[2];
sx q[2];
rz(1.3502094) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7550274) q[1];
sx q[1];
rz(-1.3530827) q[1];
sx q[1];
rz(0.45611195) q[1];
rz(2.1717908) q[3];
sx q[3];
rz(-1.9327812) q[3];
sx q[3];
rz(-2.8014744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.572523) q[2];
sx q[2];
rz(-1.5745682) q[2];
sx q[2];
rz(-1.8756867) q[2];
rz(-1.8484533) q[3];
sx q[3];
rz(-1.061941) q[3];
sx q[3];
rz(0.40614793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1311998) q[0];
sx q[0];
rz(-2.4507903) q[0];
sx q[0];
rz(-2.3660085) q[0];
rz(-0.56397437) q[1];
sx q[1];
rz(-1.0697983) q[1];
sx q[1];
rz(-1.0615798) q[1];
rz(-1.3971267) q[2];
sx q[2];
rz(-2.7026855) q[2];
sx q[2];
rz(2.4365946) q[2];
rz(2.0502144) q[3];
sx q[3];
rz(-1.8664843) q[3];
sx q[3];
rz(-2.233212) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
