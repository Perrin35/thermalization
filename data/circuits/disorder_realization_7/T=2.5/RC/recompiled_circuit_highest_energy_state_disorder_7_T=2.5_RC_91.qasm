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
rz(-1.530175) q[0];
sx q[0];
rz(2.8443008) q[0];
rz(-2.3625506) q[1];
sx q[1];
rz(-2.9809597) q[1];
sx q[1];
rz(1.7056486) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6823079) q[0];
sx q[0];
rz(-0.78792446) q[0];
sx q[0];
rz(-0.81972229) q[0];
rz(-pi) q[1];
rz(-2.024827) q[2];
sx q[2];
rz(-2.8548156) q[2];
sx q[2];
rz(1.3146626) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0432192) q[1];
sx q[1];
rz(-2.3465152) q[1];
sx q[1];
rz(2.6836392) q[1];
rz(-pi) q[2];
rz(-2.8332356) q[3];
sx q[3];
rz(-1.3801818) q[3];
sx q[3];
rz(-2.7558143) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.5982738) q[2];
sx q[2];
rz(-0.42793772) q[2];
sx q[2];
rz(-2.9738026) q[2];
rz(2.0134036) q[3];
sx q[3];
rz(-1.5690683) q[3];
sx q[3];
rz(-1.1892085) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.096864916) q[0];
sx q[0];
rz(-0.6701349) q[0];
sx q[0];
rz(2.0665533) q[0];
rz(-1.7754405) q[1];
sx q[1];
rz(-1.3864044) q[1];
sx q[1];
rz(-1.3137821) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4903298) q[0];
sx q[0];
rz(-1.8989854) q[0];
sx q[0];
rz(-1.8049678) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.73697097) q[2];
sx q[2];
rz(-2.8923375) q[2];
sx q[2];
rz(2.1971306) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.8344445) q[1];
sx q[1];
rz(-1.6475369) q[1];
sx q[1];
rz(0.84636023) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0079682) q[3];
sx q[3];
rz(-1.8989972) q[3];
sx q[3];
rz(1.4742873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.65110597) q[2];
sx q[2];
rz(-1.3765114) q[2];
sx q[2];
rz(-0.42438486) q[2];
rz(0.73244798) q[3];
sx q[3];
rz(-1.4832352) q[3];
sx q[3];
rz(1.463416) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3378147) q[0];
sx q[0];
rz(-2.4265899) q[0];
sx q[0];
rz(-2.6388229) q[0];
rz(-2.1979507) q[1];
sx q[1];
rz(-2.4991401) q[1];
sx q[1];
rz(-2.3263993) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.40476) q[0];
sx q[0];
rz(-1.6175999) q[0];
sx q[0];
rz(-2.9535049) q[0];
rz(-pi) q[1];
rz(-0.11932245) q[2];
sx q[2];
rz(-1.6762969) q[2];
sx q[2];
rz(-1.6245934) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.013061831) q[1];
sx q[1];
rz(-1.3088041) q[1];
sx q[1];
rz(-2.872741) q[1];
rz(-pi) q[2];
x q[2];
rz(0.56166537) q[3];
sx q[3];
rz(-2.1875854) q[3];
sx q[3];
rz(1.3514618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.0649071) q[2];
sx q[2];
rz(-1.4511329) q[2];
sx q[2];
rz(-2.0150851) q[2];
rz(1.7846151) q[3];
sx q[3];
rz(-2.6419736) q[3];
sx q[3];
rz(0.13993851) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.55946881) q[0];
sx q[0];
rz(-0.72002763) q[0];
sx q[0];
rz(0.35991392) q[0];
rz(0.24078807) q[1];
sx q[1];
rz(-1.1477995) q[1];
sx q[1];
rz(0.37318939) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8561962) q[0];
sx q[0];
rz(-1.4715949) q[0];
sx q[0];
rz(1.4267992) q[0];
x q[1];
rz(-0.56314205) q[2];
sx q[2];
rz(-1.3712607) q[2];
sx q[2];
rz(0.89479252) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.2149005) q[1];
sx q[1];
rz(-2.5040031) q[1];
sx q[1];
rz(-2.5465148) q[1];
x q[2];
rz(-0.4469965) q[3];
sx q[3];
rz(-1.3986386) q[3];
sx q[3];
rz(-0.76154532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.0571664) q[2];
sx q[2];
rz(-1.7095704) q[2];
sx q[2];
rz(-2.2331494) q[2];
rz(-1.069979) q[3];
sx q[3];
rz(-1.9570276) q[3];
sx q[3];
rz(-2.158304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1835566) q[0];
sx q[0];
rz(-2.2479842) q[0];
sx q[0];
rz(-2.2285484) q[0];
rz(-2.4519582) q[1];
sx q[1];
rz(-2.2328186) q[1];
sx q[1];
rz(0.77401179) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.000074) q[0];
sx q[0];
rz(-0.98709269) q[0];
sx q[0];
rz(1.0144177) q[0];
x q[1];
rz(0.795587) q[2];
sx q[2];
rz(-1.5930111) q[2];
sx q[2];
rz(-1.3327679) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.23394984) q[1];
sx q[1];
rz(-1.3798668) q[1];
sx q[1];
rz(-2.0067061) q[1];
x q[2];
rz(1.5411395) q[3];
sx q[3];
rz(-1.0231442) q[3];
sx q[3];
rz(2.1911007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.8959877) q[2];
sx q[2];
rz(-1.3934803) q[2];
sx q[2];
rz(2.4020933) q[2];
rz(1.5064404) q[3];
sx q[3];
rz(-0.93540257) q[3];
sx q[3];
rz(2.3338649) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5189811) q[0];
sx q[0];
rz(-0.7951355) q[0];
sx q[0];
rz(2.0163037) q[0];
rz(0.13825026) q[1];
sx q[1];
rz(-1.6743276) q[1];
sx q[1];
rz(1.5743871) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7049887) q[0];
sx q[0];
rz(-1.2452176) q[0];
sx q[0];
rz(2.7537936) q[0];
rz(0.41344686) q[2];
sx q[2];
rz(-1.8395429) q[2];
sx q[2];
rz(-0.15971249) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.7789845) q[1];
sx q[1];
rz(-1.4464738) q[1];
sx q[1];
rz(-2.881114) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3732143) q[3];
sx q[3];
rz(-2.1113331) q[3];
sx q[3];
rz(-0.61321248) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3681727) q[2];
sx q[2];
rz(-2.07351) q[2];
sx q[2];
rz(1.0991905) q[2];
rz(-0.98569551) q[3];
sx q[3];
rz(-1.6253977) q[3];
sx q[3];
rz(-0.43258468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2987591) q[0];
sx q[0];
rz(-2.1465813) q[0];
sx q[0];
rz(0.80879912) q[0];
rz(-2.47593) q[1];
sx q[1];
rz(-2.1620965) q[1];
sx q[1];
rz(2.1601802) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27941376) q[0];
sx q[0];
rz(-1.0004603) q[0];
sx q[0];
rz(-3.0736607) q[0];
rz(2.7000325) q[2];
sx q[2];
rz(-2.5924263) q[2];
sx q[2];
rz(1.5873945) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7715523) q[1];
sx q[1];
rz(-1.1070894) q[1];
sx q[1];
rz(0.15180219) q[1];
rz(-0.89831619) q[3];
sx q[3];
rz(-1.4631464) q[3];
sx q[3];
rz(-3.0458991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.0770646) q[2];
sx q[2];
rz(-1.4746102) q[2];
sx q[2];
rz(2.7725753) q[2];
rz(-1.7391694) q[3];
sx q[3];
rz(-2.2544421) q[3];
sx q[3];
rz(1.8203452) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8333261) q[0];
sx q[0];
rz(-0.82887355) q[0];
sx q[0];
rz(-0.33018026) q[0];
rz(-2.6847243) q[1];
sx q[1];
rz(-0.43532443) q[1];
sx q[1];
rz(0.59828573) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.24007455) q[0];
sx q[0];
rz(-0.80722729) q[0];
sx q[0];
rz(-0.62348311) q[0];
x q[1];
rz(0.070934709) q[2];
sx q[2];
rz(-2.0389551) q[2];
sx q[2];
rz(-0.93229655) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.62245) q[1];
sx q[1];
rz(-2.7323117) q[1];
sx q[1];
rz(-0.60731387) q[1];
rz(-pi) q[2];
rz(-0.71514327) q[3];
sx q[3];
rz(-0.79350457) q[3];
sx q[3];
rz(-1.8737829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2716219) q[2];
sx q[2];
rz(-2.2651256) q[2];
sx q[2];
rz(3.0769707) q[2];
rz(1.6914852) q[3];
sx q[3];
rz(-2.7414069) q[3];
sx q[3];
rz(-1.5959285) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2794063) q[0];
sx q[0];
rz(-2.8074844) q[0];
sx q[0];
rz(-2.620328) q[0];
rz(-2.0817256) q[1];
sx q[1];
rz(-1.1618549) q[1];
sx q[1];
rz(-1.0313787) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9654664) q[0];
sx q[0];
rz(-0.83399189) q[0];
sx q[0];
rz(1.470572) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6907482) q[2];
sx q[2];
rz(-1.2799096) q[2];
sx q[2];
rz(-1.6294668) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.65341144) q[1];
sx q[1];
rz(-2.0926884) q[1];
sx q[1];
rz(2.8074991) q[1];
rz(-pi) q[2];
rz(0.1886173) q[3];
sx q[3];
rz(-2.6572356) q[3];
sx q[3];
rz(0.25121197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.7445659) q[2];
sx q[2];
rz(-0.72781813) q[2];
sx q[2];
rz(3.0542206) q[2];
rz(-1.9258063) q[3];
sx q[3];
rz(-1.6814503) q[3];
sx q[3];
rz(-0.17010918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4755197) q[0];
sx q[0];
rz(-2.2505794) q[0];
sx q[0];
rz(-1.1391621) q[0];
rz(0.49925223) q[1];
sx q[1];
rz(-1.4491932) q[1];
sx q[1];
rz(-1.3334691) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.069786929) q[0];
sx q[0];
rz(-0.13955586) q[0];
sx q[0];
rz(2.4371171) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3158445) q[2];
sx q[2];
rz(-1.6189908) q[2];
sx q[2];
rz(-2.9914837) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.51133174) q[1];
sx q[1];
rz(-1.3553265) q[1];
sx q[1];
rz(0.25863077) q[1];
rz(-pi) q[2];
rz(-0.74204294) q[3];
sx q[3];
rz(-0.84845209) q[3];
sx q[3];
rz(2.6076041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.501005) q[2];
sx q[2];
rz(-1.86684) q[2];
sx q[2];
rz(-0.28323832) q[2];
rz(1.2011652) q[3];
sx q[3];
rz(-2.4656506) q[3];
sx q[3];
rz(1.0454319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.82861154) q[0];
sx q[0];
rz(-2.3647478) q[0];
sx q[0];
rz(1.1051529) q[0];
rz(-0.10440566) q[1];
sx q[1];
rz(-1.5515635) q[1];
sx q[1];
rz(-1.9952231) q[1];
rz(2.7945065) q[2];
sx q[2];
rz(-2.6095359) q[2];
sx q[2];
rz(-0.92892854) q[2];
rz(-2.6324568) q[3];
sx q[3];
rz(-0.73496277) q[3];
sx q[3];
rz(-3.1258813) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
