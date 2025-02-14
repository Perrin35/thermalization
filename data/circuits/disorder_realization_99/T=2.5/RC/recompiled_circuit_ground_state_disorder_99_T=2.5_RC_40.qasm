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
rz(-0.39781308) q[0];
sx q[0];
rz(1.8737268) q[0];
rz(-0.59029382) q[1];
sx q[1];
rz(-0.78650147) q[1];
sx q[1];
rz(-1.2193349) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7375921) q[0];
sx q[0];
rz(-0.9989748) q[0];
sx q[0];
rz(-0.19947995) q[0];
x q[1];
rz(-0.27094312) q[2];
sx q[2];
rz(-0.56519714) q[2];
sx q[2];
rz(2.924356) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.75324149) q[1];
sx q[1];
rz(-1.2920391) q[1];
sx q[1];
rz(-0.21215794) q[1];
rz(-pi) q[2];
rz(-2.9079535) q[3];
sx q[3];
rz(-0.92687449) q[3];
sx q[3];
rz(1.131191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.56403247) q[2];
sx q[2];
rz(-1.6200248) q[2];
sx q[2];
rz(-1.3684028) q[2];
rz(2.2300301) q[3];
sx q[3];
rz(-1.7715958) q[3];
sx q[3];
rz(-3.1172359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1034705) q[0];
sx q[0];
rz(-0.39762527) q[0];
sx q[0];
rz(-3.0225515) q[0];
rz(1.1301522) q[1];
sx q[1];
rz(-2.1739013) q[1];
sx q[1];
rz(1.4281323) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2909236) q[0];
sx q[0];
rz(-2.4570298) q[0];
sx q[0];
rz(-2.4638891) q[0];
x q[1];
rz(1.9910286) q[2];
sx q[2];
rz(-1.2017851) q[2];
sx q[2];
rz(-0.29847658) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.94257689) q[1];
sx q[1];
rz(-1.4369094) q[1];
sx q[1];
rz(-1.1034637) q[1];
rz(-0.30562206) q[3];
sx q[3];
rz(-1.8942617) q[3];
sx q[3];
rz(2.4572069) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1044384) q[0];
sx q[0];
rz(-2.8762682) q[0];
sx q[0];
rz(1.6695439) q[0];
rz(-0.56023359) q[1];
sx q[1];
rz(-1.7543703) q[1];
sx q[1];
rz(1.4405506) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.154523) q[0];
sx q[0];
rz(-1.3366342) q[0];
sx q[0];
rz(-2.9375141) q[0];
rz(-1.57651) q[2];
sx q[2];
rz(-2.2457321) q[2];
sx q[2];
rz(-0.99944464) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.3095197) q[1];
sx q[1];
rz(-2.4997871) q[1];
sx q[1];
rz(0.9819531) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7974924) q[3];
sx q[3];
rz(-1.3406702) q[3];
sx q[3];
rz(0.067148681) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.4855839) q[2];
sx q[2];
rz(-1.6501082) q[2];
sx q[2];
rz(-0.34240016) q[2];
rz(1.418669) q[3];
sx q[3];
rz(-2.4125621) q[3];
sx q[3];
rz(0.5595783) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2320084) q[0];
sx q[0];
rz(-0.37264687) q[0];
sx q[0];
rz(-1.6147856) q[0];
rz(-2.3341663) q[1];
sx q[1];
rz(-1.5734943) q[1];
sx q[1];
rz(-1.2015013) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.75185173) q[0];
sx q[0];
rz(-0.83503113) q[0];
sx q[0];
rz(0.40919183) q[0];
rz(2.4929659) q[2];
sx q[2];
rz(-1.0264215) q[2];
sx q[2];
rz(2.8971162) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.8462584) q[1];
sx q[1];
rz(-2.3148119) q[1];
sx q[1];
rz(2.1200402) q[1];
rz(-pi) q[2];
rz(-0.25629429) q[3];
sx q[3];
rz(-1.5206889) q[3];
sx q[3];
rz(2.9378892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3379007) q[2];
sx q[2];
rz(-0.98946977) q[2];
sx q[2];
rz(1.9039512) q[2];
rz(1.3112274) q[3];
sx q[3];
rz(-1.6236191) q[3];
sx q[3];
rz(-1.9633912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.118498) q[0];
sx q[0];
rz(-1.3730405) q[0];
sx q[0];
rz(-0.9088687) q[0];
rz(-2.2591649) q[1];
sx q[1];
rz(-2.4274223) q[1];
sx q[1];
rz(2.4095101) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1279739) q[0];
sx q[0];
rz(-0.92500702) q[0];
sx q[0];
rz(2.5396106) q[0];
rz(2.3734748) q[2];
sx q[2];
rz(-1.7752194) q[2];
sx q[2];
rz(0.26917514) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.8632224) q[1];
sx q[1];
rz(-1.178982) q[1];
sx q[1];
rz(-2.6696221) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3644712) q[3];
sx q[3];
rz(-0.91360559) q[3];
sx q[3];
rz(1.7390206) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0714134) q[2];
sx q[2];
rz(-0.81255239) q[2];
sx q[2];
rz(1.2893691) q[2];
rz(1.4687126) q[3];
sx q[3];
rz(-1.2082992) q[3];
sx q[3];
rz(-0.41951352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5656723) q[0];
sx q[0];
rz(-1.9603632) q[0];
sx q[0];
rz(2.0400203) q[0];
rz(-0.22483243) q[1];
sx q[1];
rz(-1.3460527) q[1];
sx q[1];
rz(-2.1563931) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12041872) q[0];
sx q[0];
rz(-2.2024931) q[0];
sx q[0];
rz(2.9178502) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9198138) q[2];
sx q[2];
rz(-2.1563081) q[2];
sx q[2];
rz(1.7988009) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8781153) q[1];
sx q[1];
rz(-2.2190071) q[1];
sx q[1];
rz(-0.99261673) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5959625) q[3];
sx q[3];
rz(-1.0428793) q[3];
sx q[3];
rz(2.2759469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.3605986) q[2];
sx q[2];
rz(-0.66086078) q[2];
sx q[2];
rz(-1.9471656) q[2];
rz(3.0418975) q[3];
sx q[3];
rz(-1.7710779) q[3];
sx q[3];
rz(-0.97894198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4286956) q[0];
sx q[0];
rz(-0.45658699) q[0];
sx q[0];
rz(-1.1140484) q[0];
rz(1.3975337) q[1];
sx q[1];
rz(-1.7618529) q[1];
sx q[1];
rz(2.8278415) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2148967) q[0];
sx q[0];
rz(-0.29698661) q[0];
sx q[0];
rz(-1.7420705) q[0];
x q[1];
rz(2.9767562) q[2];
sx q[2];
rz(-1.4408752) q[2];
sx q[2];
rz(2.8292953) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.3849029) q[1];
sx q[1];
rz(-2.3339565) q[1];
sx q[1];
rz(-2.7212423) q[1];
rz(0.31734316) q[3];
sx q[3];
rz(-0.71094027) q[3];
sx q[3];
rz(-1.2430354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.0745915) q[2];
sx q[2];
rz(-2.7374697) q[2];
sx q[2];
rz(-2.5355549) q[2];
rz(-2.0634985) q[3];
sx q[3];
rz(-2.2793016) q[3];
sx q[3];
rz(-1.3585453) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9692877) q[0];
sx q[0];
rz(-0.91738874) q[0];
sx q[0];
rz(-2.0148) q[0];
rz(-0.91398319) q[1];
sx q[1];
rz(-2.7028449) q[1];
sx q[1];
rz(0.21805683) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1870616) q[0];
sx q[0];
rz(-1.1370704) q[0];
sx q[0];
rz(-2.0824672) q[0];
rz(-pi) q[1];
rz(-2.281222) q[2];
sx q[2];
rz(-2.2220774) q[2];
sx q[2];
rz(1.8818579) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.14971241) q[1];
sx q[1];
rz(-2.212095) q[1];
sx q[1];
rz(-0.15705697) q[1];
rz(-1.7084684) q[3];
sx q[3];
rz(-1.6873423) q[3];
sx q[3];
rz(-2.2046409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.6286991) q[2];
sx q[2];
rz(-2.4662374) q[2];
sx q[2];
rz(2.8524354) q[2];
rz(2.1469927) q[3];
sx q[3];
rz(-1.9528439) q[3];
sx q[3];
rz(-2.6836256) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4299803) q[0];
sx q[0];
rz(-2.0531605) q[0];
sx q[0];
rz(2.3751538) q[0];
rz(1.6038766) q[1];
sx q[1];
rz(-1.8285373) q[1];
sx q[1];
rz(0.91805735) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9452451) q[0];
sx q[0];
rz(-2.1736988) q[0];
sx q[0];
rz(-2.5640998) q[0];
rz(-pi) q[1];
rz(-2.9730498) q[2];
sx q[2];
rz(-1.4282116) q[2];
sx q[2];
rz(-3.0592205) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.854771) q[1];
sx q[1];
rz(-0.76463759) q[1];
sx q[1];
rz(0.36061339) q[1];
rz(-pi) q[2];
rz(-3.0357846) q[3];
sx q[3];
rz(-1.8195767) q[3];
sx q[3];
rz(-0.63423587) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2406769) q[2];
sx q[2];
rz(-2.675481) q[2];
sx q[2];
rz(3.0898928) q[2];
rz(-2.2895571) q[3];
sx q[3];
rz(-1.7107191) q[3];
sx q[3];
rz(2.8813072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8155415) q[0];
sx q[0];
rz(-1.4651848) q[0];
sx q[0];
rz(-3.0503804) q[0];
rz(-0.7971898) q[1];
sx q[1];
rz(-1.3858831) q[1];
sx q[1];
rz(0.43112722) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1108284) q[0];
sx q[0];
rz(-1.618396) q[0];
sx q[0];
rz(3.0831227) q[0];
x q[1];
rz(1.4412142) q[2];
sx q[2];
rz(-0.65276546) q[2];
sx q[2];
rz(1.628794) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.28980892) q[1];
sx q[1];
rz(-2.0153592) q[1];
sx q[1];
rz(-1.8123933) q[1];
rz(-pi) q[2];
rz(-0.96980181) q[3];
sx q[3];
rz(-1.9327812) q[3];
sx q[3];
rz(0.34011823) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.572523) q[2];
sx q[2];
rz(-1.5670245) q[2];
sx q[2];
rz(1.265906) q[2];
rz(1.2931394) q[3];
sx q[3];
rz(-1.061941) q[3];
sx q[3];
rz(0.40614793) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1311998) q[0];
sx q[0];
rz(-0.69080234) q[0];
sx q[0];
rz(0.77558415) q[0];
rz(-2.5776183) q[1];
sx q[1];
rz(-2.0717944) q[1];
sx q[1];
rz(2.0800128) q[1];
rz(-1.7444659) q[2];
sx q[2];
rz(-0.43890719) q[2];
sx q[2];
rz(-0.70499805) q[2];
rz(-1.0913783) q[3];
sx q[3];
rz(-1.8664843) q[3];
sx q[3];
rz(-2.233212) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
