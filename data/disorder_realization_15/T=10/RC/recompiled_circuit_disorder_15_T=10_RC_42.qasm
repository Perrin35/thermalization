OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.8437334) q[0];
sx q[0];
rz(-0.61367404) q[0];
sx q[0];
rz(-2.4224129) q[0];
rz(1.367388) q[1];
sx q[1];
rz(2.8957638) q[1];
sx q[1];
rz(10.413269) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.91627097) q[0];
sx q[0];
rz(-2.7164408) q[0];
sx q[0];
rz(-1.7122373) q[0];
rz(-pi) q[1];
rz(2.9814331) q[2];
sx q[2];
rz(-0.9848435) q[2];
sx q[2];
rz(1.1908659) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.72585427) q[1];
sx q[1];
rz(-1.1032747) q[1];
sx q[1];
rz(1.0456677) q[1];
rz(-pi) q[2];
x q[2];
rz(0.037976102) q[3];
sx q[3];
rz(-1.8137323) q[3];
sx q[3];
rz(1.9248885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.44895479) q[2];
sx q[2];
rz(-1.8779034) q[2];
sx q[2];
rz(2.5732178) q[2];
rz(1.189399) q[3];
sx q[3];
rz(-0.21681924) q[3];
sx q[3];
rz(-2.9590759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi/2) q[3];
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
rz(2.4028567) q[0];
sx q[0];
rz(-1.0915272) q[0];
sx q[0];
rz(-0.36303315) q[0];
rz(2.1733213) q[1];
sx q[1];
rz(-0.67494154) q[1];
sx q[1];
rz(1.8889069) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.90447146) q[0];
sx q[0];
rz(-1.7303109) q[0];
sx q[0];
rz(-1.479854) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.91141869) q[2];
sx q[2];
rz(-1.0936001) q[2];
sx q[2];
rz(-2.0267682) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.82488793) q[1];
sx q[1];
rz(-2.0522293) q[1];
sx q[1];
rz(0.96932051) q[1];
rz(1.1329123) q[3];
sx q[3];
rz(-1.7403733) q[3];
sx q[3];
rz(1.358043) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.12053717) q[2];
sx q[2];
rz(-1.8214104) q[2];
sx q[2];
rz(-2.7978314) q[2];
rz(2.8072642) q[3];
sx q[3];
rz(-0.40083405) q[3];
sx q[3];
rz(-0.46935579) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47135982) q[0];
sx q[0];
rz(-2.8911262) q[0];
sx q[0];
rz(0.0070455889) q[0];
rz(0.37653157) q[1];
sx q[1];
rz(-2.2129602) q[1];
sx q[1];
rz(0.71281707) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.26582742) q[0];
sx q[0];
rz(-1.3607425) q[0];
sx q[0];
rz(1.34583) q[0];
rz(-1.1586645) q[2];
sx q[2];
rz(-2.7325028) q[2];
sx q[2];
rz(-1.2661753) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.8411326) q[1];
sx q[1];
rz(-1.7736048) q[1];
sx q[1];
rz(2.1122785) q[1];
rz(-pi) q[2];
rz(-2.0531822) q[3];
sx q[3];
rz(-1.6420806) q[3];
sx q[3];
rz(0.34382581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.787848) q[2];
sx q[2];
rz(-1.5622082) q[2];
sx q[2];
rz(0.66398579) q[2];
rz(1.9021696) q[3];
sx q[3];
rz(-2.911471) q[3];
sx q[3];
rz(-0.91528875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5749213) q[0];
sx q[0];
rz(-3.1118588) q[0];
sx q[0];
rz(0.57408875) q[0];
rz(-0.29218778) q[1];
sx q[1];
rz(-0.27583396) q[1];
sx q[1];
rz(-0.92837292) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.95604491) q[0];
sx q[0];
rz(-1.569824) q[0];
sx q[0];
rz(2.2768343) q[0];
rz(-pi) q[1];
x q[1];
rz(0.78903918) q[2];
sx q[2];
rz(-1.3256729) q[2];
sx q[2];
rz(2.6546728) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1410364) q[1];
sx q[1];
rz(-2.6489241) q[1];
sx q[1];
rz(2.1129184) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2286962) q[3];
sx q[3];
rz(-1.9808589) q[3];
sx q[3];
rz(2.0994772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.092768) q[2];
sx q[2];
rz(-0.51468819) q[2];
sx q[2];
rz(-1.9862004) q[2];
rz(-2.998735) q[3];
sx q[3];
rz(-2.6070049) q[3];
sx q[3];
rz(0.75240451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.083754152) q[0];
sx q[0];
rz(-2.6733119) q[0];
sx q[0];
rz(-3.0673448) q[0];
rz(1.4986562) q[1];
sx q[1];
rz(-1.5131806) q[1];
sx q[1];
rz(0.9517076) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4518406) q[0];
sx q[0];
rz(-1.6432439) q[0];
sx q[0];
rz(-0.82093443) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.81163002) q[2];
sx q[2];
rz(-1.1434165) q[2];
sx q[2];
rz(0.4610093) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.4068027) q[1];
sx q[1];
rz(-2.4996335) q[1];
sx q[1];
rz(1.2924679) q[1];
rz(-pi) q[2];
rz(-2.1390901) q[3];
sx q[3];
rz(-1.2270524) q[3];
sx q[3];
rz(-0.4337726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4962861) q[2];
sx q[2];
rz(-0.22660613) q[2];
sx q[2];
rz(-1.9892233) q[2];
rz(1.0460098) q[3];
sx q[3];
rz(-0.73863107) q[3];
sx q[3];
rz(-0.0013105198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1482658) q[0];
sx q[0];
rz(-0.95411211) q[0];
sx q[0];
rz(2.6519725) q[0];
rz(2.1173677) q[1];
sx q[1];
rz(-2.5074597) q[1];
sx q[1];
rz(1.6061868) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0438457) q[0];
sx q[0];
rz(-1.6698208) q[0];
sx q[0];
rz(1.6188341) q[0];
rz(-2.4602731) q[2];
sx q[2];
rz(-1.0208924) q[2];
sx q[2];
rz(0.51598179) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.225243) q[1];
sx q[1];
rz(-0.32685977) q[1];
sx q[1];
rz(2.0845695) q[1];
x q[2];
rz(-0.225004) q[3];
sx q[3];
rz(-0.64151728) q[3];
sx q[3];
rz(1.0596421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.17807047) q[2];
sx q[2];
rz(-2.3242951) q[2];
sx q[2];
rz(2.9658588) q[2];
rz(1.0774311) q[3];
sx q[3];
rz(-1.146799) q[3];
sx q[3];
rz(3.1350465) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.068129383) q[0];
sx q[0];
rz(-0.21502762) q[0];
sx q[0];
rz(1.2699132) q[0];
rz(-2.9636256) q[1];
sx q[1];
rz(-2.3033419) q[1];
sx q[1];
rz(1.7756745) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38458347) q[0];
sx q[0];
rz(-0.88338822) q[0];
sx q[0];
rz(-1.13152) q[0];
rz(-pi) q[1];
rz(1.3283417) q[2];
sx q[2];
rz(-1.2521267) q[2];
sx q[2];
rz(0.53675011) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.3733702) q[1];
sx q[1];
rz(-1.7473514) q[1];
sx q[1];
rz(2.8819487) q[1];
x q[2];
rz(1.4468206) q[3];
sx q[3];
rz(-1.0448536) q[3];
sx q[3];
rz(-0.40243173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0718677) q[2];
sx q[2];
rz(-2.8475927) q[2];
sx q[2];
rz(-0.89795566) q[2];
rz(-0.14363025) q[3];
sx q[3];
rz(-1.2575905) q[3];
sx q[3];
rz(2.7974421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1865858) q[0];
sx q[0];
rz(-0.85383677) q[0];
sx q[0];
rz(2.8198077) q[0];
rz(0.92542648) q[1];
sx q[1];
rz(-1.7089475) q[1];
sx q[1];
rz(0.61703533) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2935534) q[0];
sx q[0];
rz(-2.8176421) q[0];
sx q[0];
rz(-0.89963161) q[0];
rz(-1.4204942) q[2];
sx q[2];
rz(-0.21810025) q[2];
sx q[2];
rz(0.66767603) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.9279234) q[1];
sx q[1];
rz(-0.59487768) q[1];
sx q[1];
rz(0.87022123) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5259605) q[3];
sx q[3];
rz(-2.3196844) q[3];
sx q[3];
rz(0.84079784) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.59163219) q[2];
sx q[2];
rz(-0.98667163) q[2];
sx q[2];
rz(-0.11403306) q[2];
rz(-2.7791801) q[3];
sx q[3];
rz(-0.24505469) q[3];
sx q[3];
rz(-1.2362278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6645901) q[0];
sx q[0];
rz(-2.6053071) q[0];
sx q[0];
rz(-1.3569008) q[0];
rz(0.82018954) q[1];
sx q[1];
rz(-0.35239041) q[1];
sx q[1];
rz(1.483451) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9634919) q[0];
sx q[0];
rz(-0.083847001) q[0];
sx q[0];
rz(1.6944147) q[0];
x q[1];
rz(-1.4610897) q[2];
sx q[2];
rz(-0.98329558) q[2];
sx q[2];
rz(2.7320931) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.2918313) q[1];
sx q[1];
rz(-1.2408248) q[1];
sx q[1];
rz(-2.989819) q[1];
x q[2];
rz(-1.9547192) q[3];
sx q[3];
rz(-2.0878289) q[3];
sx q[3];
rz(-2.8843055) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0004398) q[2];
sx q[2];
rz(-2.4856462) q[2];
sx q[2];
rz(1.2443939) q[2];
rz(-0.42516285) q[3];
sx q[3];
rz(-0.77912283) q[3];
sx q[3];
rz(-1.6311133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4193831) q[0];
sx q[0];
rz(-0.49383759) q[0];
sx q[0];
rz(2.7710932) q[0];
rz(-2.0314979) q[1];
sx q[1];
rz(-0.67567527) q[1];
sx q[1];
rz(1.8006905) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6275218) q[0];
sx q[0];
rz(-0.43826575) q[0];
sx q[0];
rz(-2.1379495) q[0];
rz(-pi) q[1];
rz(2.9843763) q[2];
sx q[2];
rz(-0.72050691) q[2];
sx q[2];
rz(2.9027028) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1738051) q[1];
sx q[1];
rz(-0.44631821) q[1];
sx q[1];
rz(2.1800969) q[1];
x q[2];
rz(3.020346) q[3];
sx q[3];
rz(-1.5505425) q[3];
sx q[3];
rz(-2.1840087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.5122539) q[2];
sx q[2];
rz(-0.24015716) q[2];
sx q[2];
rz(1.5226927) q[2];
rz(2.2807138) q[3];
sx q[3];
rz(-1.7247793) q[3];
sx q[3];
rz(0.035877429) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2387977) q[0];
sx q[0];
rz(-1.6607941) q[0];
sx q[0];
rz(-0.86097062) q[0];
rz(0.32342708) q[1];
sx q[1];
rz(-1.885781) q[1];
sx q[1];
rz(1.5992004) q[1];
rz(1.1644438) q[2];
sx q[2];
rz(-1.5487557) q[2];
sx q[2];
rz(1.1974481) q[2];
rz(-2.042751) q[3];
sx q[3];
rz(-1.9190211) q[3];
sx q[3];
rz(2.2091051) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
