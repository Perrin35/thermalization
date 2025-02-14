OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.4647291) q[0];
sx q[0];
rz(1.6114177) q[0];
sx q[0];
rz(9.7220698) q[0];
rz(-2.3625506) q[1];
sx q[1];
rz(-2.9809597) q[1];
sx q[1];
rz(1.7056486) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6061062) q[0];
sx q[0];
rz(-2.1155042) q[0];
sx q[0];
rz(-0.60118712) q[0];
x q[1];
rz(2.024827) q[2];
sx q[2];
rz(-2.8548156) q[2];
sx q[2];
rz(-1.3146626) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.7117225) q[1];
sx q[1];
rz(-2.2657569) q[1];
sx q[1];
rz(1.9942787) q[1];
x q[2];
rz(-1.3709896) q[3];
sx q[3];
rz(-1.268203) q[3];
sx q[3];
rz(-1.8963008) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5982738) q[2];
sx q[2];
rz(-2.7136549) q[2];
sx q[2];
rz(0.16779009) q[2];
rz(-2.0134036) q[3];
sx q[3];
rz(-1.5690683) q[3];
sx q[3];
rz(1.1892085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0447277) q[0];
sx q[0];
rz(-0.6701349) q[0];
sx q[0];
rz(1.0750394) q[0];
rz(1.3661522) q[1];
sx q[1];
rz(-1.3864044) q[1];
sx q[1];
rz(1.8278106) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2883817) q[0];
sx q[0];
rz(-0.40069212) q[0];
sx q[0];
rz(0.59817065) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.4013675) q[2];
sx q[2];
rz(-1.3870948) q[2];
sx q[2];
rz(1.6970762) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.35006501) q[1];
sx q[1];
rz(-0.72775562) q[1];
sx q[1];
rz(-1.6863053) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9437378) q[3];
sx q[3];
rz(-2.7881456) q[3];
sx q[3];
rz(1.8693876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.4904867) q[2];
sx q[2];
rz(-1.7650812) q[2];
sx q[2];
rz(2.7172078) q[2];
rz(-2.4091447) q[3];
sx q[3];
rz(-1.6583574) q[3];
sx q[3];
rz(1.6781767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80377793) q[0];
sx q[0];
rz(-0.71500272) q[0];
sx q[0];
rz(-2.6388229) q[0];
rz(0.94364199) q[1];
sx q[1];
rz(-0.6424526) q[1];
sx q[1];
rz(2.3263993) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0749933) q[0];
sx q[0];
rz(-0.19375676) q[0];
sx q[0];
rz(-2.8961477) q[0];
rz(-0.72702144) q[2];
sx q[2];
rz(-2.9824848) q[2];
sx q[2];
rz(-2.4746759) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(3.1285308) q[1];
sx q[1];
rz(-1.3088041) q[1];
sx q[1];
rz(2.872741) q[1];
rz(-pi) q[2];
rz(2.5799273) q[3];
sx q[3];
rz(-2.1875854) q[3];
sx q[3];
rz(1.7901309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.076685585) q[2];
sx q[2];
rz(-1.4511329) q[2];
sx q[2];
rz(1.1265075) q[2];
rz(-1.7846151) q[3];
sx q[3];
rz(-2.6419736) q[3];
sx q[3];
rz(-0.13993851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55946881) q[0];
sx q[0];
rz(-2.421565) q[0];
sx q[0];
rz(2.7816787) q[0];
rz(0.24078807) q[1];
sx q[1];
rz(-1.9937932) q[1];
sx q[1];
rz(-0.37318939) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2853965) q[0];
sx q[0];
rz(-1.4715949) q[0];
sx q[0];
rz(-1.7147934) q[0];
x q[1];
rz(2.5784506) q[2];
sx q[2];
rz(-1.7703319) q[2];
sx q[2];
rz(-0.89479252) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.2265985) q[1];
sx q[1];
rz(-1.0553331) q[1];
sx q[1];
rz(1.9644009) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.38244989) q[3];
sx q[3];
rz(-2.6646864) q[3];
sx q[3];
rz(-0.46602369) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0571664) q[2];
sx q[2];
rz(-1.7095704) q[2];
sx q[2];
rz(-2.2331494) q[2];
rz(1.069979) q[3];
sx q[3];
rz(-1.9570276) q[3];
sx q[3];
rz(2.158304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95803607) q[0];
sx q[0];
rz(-0.89360845) q[0];
sx q[0];
rz(0.91304427) q[0];
rz(2.4519582) q[1];
sx q[1];
rz(-2.2328186) q[1];
sx q[1];
rz(2.3675809) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1415187) q[0];
sx q[0];
rz(-0.98709269) q[0];
sx q[0];
rz(2.1271749) q[0];
rz(-pi) q[1];
rz(-0.795587) q[2];
sx q[2];
rz(-1.5930111) q[2];
sx q[2];
rz(-1.8088248) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.2486813) q[1];
sx q[1];
rz(-1.9982583) q[1];
sx q[1];
rz(0.21007725) q[1];
x q[2];
rz(2.5937449) q[3];
sx q[3];
rz(-1.5454779) q[3];
sx q[3];
rz(2.5058432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.245605) q[2];
sx q[2];
rz(-1.3934803) q[2];
sx q[2];
rz(2.4020933) q[2];
rz(1.6351522) q[3];
sx q[3];
rz(-0.93540257) q[3];
sx q[3];
rz(-2.3338649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6226115) q[0];
sx q[0];
rz(-2.3464572) q[0];
sx q[0];
rz(-1.125289) q[0];
rz(0.13825026) q[1];
sx q[1];
rz(-1.467265) q[1];
sx q[1];
rz(-1.5743871) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.137321) q[0];
sx q[0];
rz(-1.204365) q[0];
sx q[0];
rz(-1.9204813) q[0];
rz(-pi) q[1];
rz(-0.41344686) q[2];
sx q[2];
rz(-1.3020497) q[2];
sx q[2];
rz(2.9818802) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.36260819) q[1];
sx q[1];
rz(-1.6951188) q[1];
sx q[1];
rz(-2.881114) q[1];
rz(-0.31611021) q[3];
sx q[3];
rz(-0.5721285) q[3];
sx q[3];
rz(-2.8993949) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.77341998) q[2];
sx q[2];
rz(-1.0680826) q[2];
sx q[2];
rz(-1.0991905) q[2];
rz(0.98569551) q[3];
sx q[3];
rz(-1.516195) q[3];
sx q[3];
rz(-0.43258468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8428335) q[0];
sx q[0];
rz(-0.99501139) q[0];
sx q[0];
rz(-2.3327935) q[0];
rz(2.47593) q[1];
sx q[1];
rz(-0.97949615) q[1];
sx q[1];
rz(2.1601802) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8869276) q[0];
sx q[0];
rz(-1.5136295) q[0];
sx q[0];
rz(0.99941038) q[0];
rz(-pi) q[1];
rz(-2.7000325) q[2];
sx q[2];
rz(-2.5924263) q[2];
sx q[2];
rz(-1.5873945) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.0404824) q[1];
sx q[1];
rz(-2.6553934) q[1];
sx q[1];
rz(-1.8644488) q[1];
x q[2];
rz(3.0043169) q[3];
sx q[3];
rz(-0.90291728) q[3];
sx q[3];
rz(1.5811416) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.064528) q[2];
sx q[2];
rz(-1.4746102) q[2];
sx q[2];
rz(-0.36901739) q[2];
rz(1.4024233) q[3];
sx q[3];
rz(-2.2544421) q[3];
sx q[3];
rz(1.8203452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3082665) q[0];
sx q[0];
rz(-2.3127191) q[0];
sx q[0];
rz(-0.33018026) q[0];
rz(0.45686832) q[1];
sx q[1];
rz(-0.43532443) q[1];
sx q[1];
rz(0.59828573) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5766524) q[0];
sx q[0];
rz(-2.1974753) q[0];
sx q[0];
rz(-2.1184854) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0706579) q[2];
sx q[2];
rz(-2.0389551) q[2];
sx q[2];
rz(2.2092961) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6574994) q[1];
sx q[1];
rz(-1.3417021) q[1];
sx q[1];
rz(-2.799396) q[1];
rz(-0.65450683) q[3];
sx q[3];
rz(-1.084436) q[3];
sx q[3];
rz(0.24412046) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2716219) q[2];
sx q[2];
rz(-2.2651256) q[2];
sx q[2];
rz(-0.064621933) q[2];
rz(-1.4501075) q[3];
sx q[3];
rz(-0.40018574) q[3];
sx q[3];
rz(-1.5456642) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2794063) q[0];
sx q[0];
rz(-2.8074844) q[0];
sx q[0];
rz(0.52126467) q[0];
rz(1.0598671) q[1];
sx q[1];
rz(-1.1618549) q[1];
sx q[1];
rz(-1.0313787) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3272031) q[0];
sx q[0];
rz(-1.6449682) q[0];
sx q[0];
rz(2.402284) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6907482) q[2];
sx q[2];
rz(-1.2799096) q[2];
sx q[2];
rz(1.6294668) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.045193521) q[1];
sx q[1];
rz(-0.61122433) q[1];
sx q[1];
rz(-1.0525714) q[1];
x q[2];
rz(2.9529753) q[3];
sx q[3];
rz(-2.6572356) q[3];
sx q[3];
rz(-0.25121197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7445659) q[2];
sx q[2];
rz(-0.72781813) q[2];
sx q[2];
rz(0.087372027) q[2];
rz(1.9258063) q[3];
sx q[3];
rz(-1.4601424) q[3];
sx q[3];
rz(2.9714835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4755197) q[0];
sx q[0];
rz(-2.2505794) q[0];
sx q[0];
rz(-2.0024306) q[0];
rz(-2.6423404) q[1];
sx q[1];
rz(-1.6923994) q[1];
sx q[1];
rz(-1.8081236) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0718057) q[0];
sx q[0];
rz(-3.0020368) q[0];
sx q[0];
rz(0.7044756) q[0];
rz(0.049801783) q[2];
sx q[2];
rz(-1.3161471) q[2];
sx q[2];
rz(-1.7083502) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.6302609) q[1];
sx q[1];
rz(-1.3553265) q[1];
sx q[1];
rz(-0.25863077) q[1];
x q[2];
rz(2.3995497) q[3];
sx q[3];
rz(-0.84845209) q[3];
sx q[3];
rz(2.6076041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.501005) q[2];
sx q[2];
rz(-1.2747526) q[2];
sx q[2];
rz(-2.8583543) q[2];
rz(-1.9404274) q[3];
sx q[3];
rz(-2.4656506) q[3];
sx q[3];
rz(-2.0961608) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3129811) q[0];
sx q[0];
rz(-0.77684488) q[0];
sx q[0];
rz(-2.0364398) q[0];
rz(3.037187) q[1];
sx q[1];
rz(-1.5515635) q[1];
sx q[1];
rz(-1.9952231) q[1];
rz(1.3731643) q[2];
sx q[2];
rz(-1.0734954) q[2];
sx q[2];
rz(-1.3263477) q[2];
rz(2.4734409) q[3];
sx q[3];
rz(-1.9037608) q[3];
sx q[3];
rz(-1.9477061) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
