OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.428838) q[0];
sx q[0];
rz(-1.5313671) q[0];
sx q[0];
rz(-1.1046326) q[0];
rz(-1.8072577) q[1];
sx q[1];
rz(-0.72307888) q[1];
sx q[1];
rz(-1.877797) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9845147) q[0];
sx q[0];
rz(-0.70467585) q[0];
sx q[0];
rz(-1.1873755) q[0];
x q[1];
rz(2.14416) q[2];
sx q[2];
rz(-1.9321529) q[2];
sx q[2];
rz(-0.45623764) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.3448779) q[1];
sx q[1];
rz(-1.895322) q[1];
sx q[1];
rz(3.0561563) q[1];
x q[2];
rz(0.24407152) q[3];
sx q[3];
rz(-1.9979949) q[3];
sx q[3];
rz(-2.0352767) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9870712) q[2];
sx q[2];
rz(-2.0646586) q[2];
sx q[2];
rz(-1.7729574) q[2];
rz(0.23855071) q[3];
sx q[3];
rz(-0.41540256) q[3];
sx q[3];
rz(0.11428782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.7411165) q[0];
sx q[0];
rz(-1.1392765) q[0];
sx q[0];
rz(1.1503295) q[0];
rz(-0.087609619) q[1];
sx q[1];
rz(-1.3299512) q[1];
sx q[1];
rz(1.5706583) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99518665) q[0];
sx q[0];
rz(-2.0530465) q[0];
sx q[0];
rz(-1.6907808) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9628062) q[2];
sx q[2];
rz(-1.6331722) q[2];
sx q[2];
rz(-1.5449926) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1024061) q[1];
sx q[1];
rz(-2.3433609) q[1];
sx q[1];
rz(2.185553) q[1];
x q[2];
rz(-2.3493166) q[3];
sx q[3];
rz(-1.918692) q[3];
sx q[3];
rz(-2.874305) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.34438434) q[2];
sx q[2];
rz(-2.4281561) q[2];
sx q[2];
rz(0.1864645) q[2];
rz(-0.79948419) q[3];
sx q[3];
rz(-1.5461642) q[3];
sx q[3];
rz(-0.98108393) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.8568273) q[0];
sx q[0];
rz(-1.7600049) q[0];
sx q[0];
rz(3.0322266) q[0];
rz(0.60375396) q[1];
sx q[1];
rz(-1.8021288) q[1];
sx q[1];
rz(-2.107479) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9700978) q[0];
sx q[0];
rz(-1.7037647) q[0];
sx q[0];
rz(1.7823969) q[0];
rz(1.500152) q[2];
sx q[2];
rz(-1.4817258) q[2];
sx q[2];
rz(-1.3180817) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.32419606) q[1];
sx q[1];
rz(-1.4651555) q[1];
sx q[1];
rz(-0.28075851) q[1];
x q[2];
rz(-2.3699371) q[3];
sx q[3];
rz(-1.4691741) q[3];
sx q[3];
rz(1.1751428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.6094531) q[2];
sx q[2];
rz(-2.1707462) q[2];
sx q[2];
rz(-2.5965221) q[2];
rz(2.3405781) q[3];
sx q[3];
rz(-2.2478734) q[3];
sx q[3];
rz(-3.0494704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48701778) q[0];
sx q[0];
rz(-1.4734522) q[0];
sx q[0];
rz(0.45904485) q[0];
rz(1.1162988) q[1];
sx q[1];
rz(-2.3479159) q[1];
sx q[1];
rz(-2.3150516) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5640663) q[0];
sx q[0];
rz(-2.5788529) q[0];
sx q[0];
rz(-0.86842815) q[0];
x q[1];
rz(2.3514296) q[2];
sx q[2];
rz(-0.1383257) q[2];
sx q[2];
rz(-0.1639072) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.7956808) q[1];
sx q[1];
rz(-1.7670146) q[1];
sx q[1];
rz(-1.0190585) q[1];
rz(-2.8119254) q[3];
sx q[3];
rz(-2.7894944) q[3];
sx q[3];
rz(0.47131495) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.7894342) q[2];
sx q[2];
rz(-0.23098478) q[2];
sx q[2];
rz(-2.3918772) q[2];
rz(-2.0189144) q[3];
sx q[3];
rz(-1.9176982) q[3];
sx q[3];
rz(0.96021715) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.66626755) q[0];
sx q[0];
rz(-0.98639494) q[0];
sx q[0];
rz(-0.112003) q[0];
rz(-0.22459596) q[1];
sx q[1];
rz(-1.9738395) q[1];
sx q[1];
rz(-0.78132838) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.76666895) q[0];
sx q[0];
rz(-1.5974853) q[0];
sx q[0];
rz(-1.1666537) q[0];
rz(-2.8113643) q[2];
sx q[2];
rz(-1.4817186) q[2];
sx q[2];
rz(1.1116127) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.71415662) q[1];
sx q[1];
rz(-0.86227741) q[1];
sx q[1];
rz(2.5664213) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1634401) q[3];
sx q[3];
rz(-1.4362122) q[3];
sx q[3];
rz(0.59813598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.8319228) q[2];
sx q[2];
rz(-2.2396542) q[2];
sx q[2];
rz(2.208948) q[2];
rz(-2.0617088) q[3];
sx q[3];
rz(-2.6974758) q[3];
sx q[3];
rz(-0.035695765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2340853) q[0];
sx q[0];
rz(-2.0103173) q[0];
sx q[0];
rz(1.3908516) q[0];
rz(2.8098409) q[1];
sx q[1];
rz(-1.876588) q[1];
sx q[1];
rz(-1.5914241) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46813477) q[0];
sx q[0];
rz(-2.4979257) q[0];
sx q[0];
rz(-1.629384) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6076902) q[2];
sx q[2];
rz(-1.1114745) q[2];
sx q[2];
rz(1.5609891) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.5065803) q[1];
sx q[1];
rz(-1.4476579) q[1];
sx q[1];
rz(-1.413373) q[1];
rz(1.9523296) q[3];
sx q[3];
rz(-1.23151) q[3];
sx q[3];
rz(1.8041704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.8288237) q[2];
sx q[2];
rz(-2.2569816) q[2];
sx q[2];
rz(-1.5213607) q[2];
rz(0.96794266) q[3];
sx q[3];
rz(-1.5827551) q[3];
sx q[3];
rz(-0.91606417) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62436002) q[0];
sx q[0];
rz(-2.7150798) q[0];
sx q[0];
rz(2.970001) q[0];
rz(1.0393556) q[1];
sx q[1];
rz(-1.8828705) q[1];
sx q[1];
rz(1.7003869) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.489451) q[0];
sx q[0];
rz(-1.891937) q[0];
sx q[0];
rz(1.4399066) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6827379) q[2];
sx q[2];
rz(-1.0575235) q[2];
sx q[2];
rz(1.7626732) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6771779) q[1];
sx q[1];
rz(-2.5485793) q[1];
sx q[1];
rz(-0.89927425) q[1];
rz(1.0073184) q[3];
sx q[3];
rz(-1.451056) q[3];
sx q[3];
rz(3.0547676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.27199304) q[2];
sx q[2];
rz(-1.7337493) q[2];
sx q[2];
rz(-2.5970411) q[2];
rz(0.49992391) q[3];
sx q[3];
rz(-2.1130424) q[3];
sx q[3];
rz(-2.669529) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2917824) q[0];
sx q[0];
rz(-0.53684679) q[0];
sx q[0];
rz(2.612402) q[0];
rz(-0.57580194) q[1];
sx q[1];
rz(-1.2628097) q[1];
sx q[1];
rz(1.1189438) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7125583) q[0];
sx q[0];
rz(-3.0858485) q[0];
sx q[0];
rz(2.7554465) q[0];
x q[1];
rz(-2.7706657) q[2];
sx q[2];
rz(-1.7417522) q[2];
sx q[2];
rz(-1.2990739) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.8994979) q[1];
sx q[1];
rz(-2.3589954) q[1];
sx q[1];
rz(2.461754) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.921245) q[3];
sx q[3];
rz(-2.5867037) q[3];
sx q[3];
rz(0.12400907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.90073663) q[2];
sx q[2];
rz(-2.6764328) q[2];
sx q[2];
rz(0.71211234) q[2];
rz(1.5322878) q[3];
sx q[3];
rz(-2.1963547) q[3];
sx q[3];
rz(1.7942662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79494548) q[0];
sx q[0];
rz(-1.1279339) q[0];
sx q[0];
rz(-1.0711063) q[0];
rz(1.2062997) q[1];
sx q[1];
rz(-2.1251528) q[1];
sx q[1];
rz(-1.6546904) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3275571) q[0];
sx q[0];
rz(-2.2326755) q[0];
sx q[0];
rz(-1.8083841) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.73545154) q[2];
sx q[2];
rz(-1.1142154) q[2];
sx q[2];
rz(0.18010715) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.35569977) q[1];
sx q[1];
rz(-1.744387) q[1];
sx q[1];
rz(-2.8462571) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8603691) q[3];
sx q[3];
rz(-2.0042002) q[3];
sx q[3];
rz(3.0126257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.6335166) q[2];
sx q[2];
rz(-0.8997007) q[2];
sx q[2];
rz(-0.077795204) q[2];
rz(-0.22732321) q[3];
sx q[3];
rz(-1.1319755) q[3];
sx q[3];
rz(-1.0556861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28854293) q[0];
sx q[0];
rz(-2.396614) q[0];
sx q[0];
rz(2.7689834) q[0];
rz(-2.695072) q[1];
sx q[1];
rz(-1.4563072) q[1];
sx q[1];
rz(2.2850697) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73222173) q[0];
sx q[0];
rz(-1.5268832) q[0];
sx q[0];
rz(0.056890566) q[0];
x q[1];
rz(0.90214731) q[2];
sx q[2];
rz(-2.9716316) q[2];
sx q[2];
rz(0.68357498) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.0622934) q[1];
sx q[1];
rz(-0.41626272) q[1];
sx q[1];
rz(-1.2757311) q[1];
x q[2];
rz(-2.2482199) q[3];
sx q[3];
rz(-0.89266289) q[3];
sx q[3];
rz(-0.03420364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.38241688) q[2];
sx q[2];
rz(-2.0903812) q[2];
sx q[2];
rz(-2.5860533) q[2];
rz(0.38604745) q[3];
sx q[3];
rz(-1.4945533) q[3];
sx q[3];
rz(-0.66421318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0059218) q[0];
sx q[0];
rz(-1.4501403) q[0];
sx q[0];
rz(-1.8474664) q[0];
rz(2.338943) q[1];
sx q[1];
rz(-1.6812656) q[1];
sx q[1];
rz(0.38801286) q[1];
rz(2.8254208) q[2];
sx q[2];
rz(-2.7718622) q[2];
sx q[2];
rz(-2.3864701) q[2];
rz(1.4925719) q[3];
sx q[3];
rz(-2.409392) q[3];
sx q[3];
rz(0.87421855) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
