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
rz(-2.7231554) q[0];
sx q[0];
rz(-2.1783481) q[0];
sx q[0];
rz(2.9377687) q[0];
rz(1.0526429) q[1];
sx q[1];
rz(3.9349603) q[1];
sx q[1];
rz(10.959672) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14102916) q[0];
sx q[0];
rz(-1.3357497) q[0];
sx q[0];
rz(-2.0268593) q[0];
x q[1];
rz(0.84759961) q[2];
sx q[2];
rz(-1.3061451) q[2];
sx q[2];
rz(-0.48993313) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.71434778) q[1];
sx q[1];
rz(-2.3934919) q[1];
sx q[1];
rz(2.1022878) q[1];
rz(-pi) q[2];
rz(0.58310572) q[3];
sx q[3];
rz(-2.8994377) q[3];
sx q[3];
rz(-0.21447578) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7242929) q[2];
sx q[2];
rz(-0.55880004) q[2];
sx q[2];
rz(-2.144045) q[2];
rz(2.3857462) q[3];
sx q[3];
rz(-1.3563124) q[3];
sx q[3];
rz(2.1006179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34870979) q[0];
sx q[0];
rz(-3.0443158) q[0];
sx q[0];
rz(1.8541699) q[0];
rz(-2.6584794) q[1];
sx q[1];
rz(-2.099642) q[1];
sx q[1];
rz(-1.9226673) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3588803) q[0];
sx q[0];
rz(-0.25280372) q[0];
sx q[0];
rz(-0.88835277) q[0];
rz(1.7886247) q[2];
sx q[2];
rz(-0.91060293) q[2];
sx q[2];
rz(-0.24581395) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7658577) q[1];
sx q[1];
rz(-1.8490144) q[1];
sx q[1];
rz(-0.15032676) q[1];
x q[2];
rz(-1.990149) q[3];
sx q[3];
rz(-2.4046728) q[3];
sx q[3];
rz(1.1518971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.37811849) q[2];
sx q[2];
rz(-3.0749622) q[2];
sx q[2];
rz(0.62180579) q[2];
rz(0.63831896) q[3];
sx q[3];
rz(-2.392605) q[3];
sx q[3];
rz(-0.84426713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3062375) q[0];
sx q[0];
rz(-2.4995646) q[0];
sx q[0];
rz(2.9252692) q[0];
rz(2.5430039) q[1];
sx q[1];
rz(-1.3052156) q[1];
sx q[1];
rz(0.52282202) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4958333) q[0];
sx q[0];
rz(-1.8765939) q[0];
sx q[0];
rz(1.9093139) q[0];
x q[1];
rz(2.7504146) q[2];
sx q[2];
rz(-1.7157946) q[2];
sx q[2];
rz(2.7716293) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.4292517) q[1];
sx q[1];
rz(-2.0852226) q[1];
sx q[1];
rz(1.8271258) q[1];
rz(-2.834972) q[3];
sx q[3];
rz(-1.2546179) q[3];
sx q[3];
rz(-0.090360377) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.9598976) q[2];
sx q[2];
rz(-1.2744224) q[2];
sx q[2];
rz(-1.6240906) q[2];
rz(-2.6767139) q[3];
sx q[3];
rz(-0.93022323) q[3];
sx q[3];
rz(1.6200861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9972123) q[0];
sx q[0];
rz(-2.7535487) q[0];
sx q[0];
rz(-1.5390747) q[0];
rz(-2.1082361) q[1];
sx q[1];
rz(-2.3732503) q[1];
sx q[1];
rz(-1.296952) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.813445) q[0];
sx q[0];
rz(-1.362365) q[0];
sx q[0];
rz(-2.3013902) q[0];
rz(3.0791666) q[2];
sx q[2];
rz(-2.2656144) q[2];
sx q[2];
rz(0.42499396) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.9819543) q[1];
sx q[1];
rz(-0.57483638) q[1];
sx q[1];
rz(-2.3153618) q[1];
rz(0.4850895) q[3];
sx q[3];
rz(-2.5311433) q[3];
sx q[3];
rz(1.4097253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3492655) q[2];
sx q[2];
rz(-2.5040864) q[2];
sx q[2];
rz(2.0959334) q[2];
rz(2.9271434) q[3];
sx q[3];
rz(-2.4172473) q[3];
sx q[3];
rz(1.1469871) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8056718) q[0];
sx q[0];
rz(-1.2821953) q[0];
sx q[0];
rz(-2.6222498) q[0];
rz(0.94201159) q[1];
sx q[1];
rz(-0.28128925) q[1];
sx q[1];
rz(1.6090144) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9056775) q[0];
sx q[0];
rz(-2.2122243) q[0];
sx q[0];
rz(2.0302613) q[0];
rz(-pi) q[1];
rz(1.7132883) q[2];
sx q[2];
rz(-1.7497369) q[2];
sx q[2];
rz(-1.6662836) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.11621257) q[1];
sx q[1];
rz(-1.0702225) q[1];
sx q[1];
rz(3.1356372) q[1];
rz(3.0318465) q[3];
sx q[3];
rz(-1.6648431) q[3];
sx q[3];
rz(1.7141683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.64451009) q[2];
sx q[2];
rz(-2.3641059) q[2];
sx q[2];
rz(2.8572148) q[2];
rz(-0.55822462) q[3];
sx q[3];
rz(-1.3642718) q[3];
sx q[3];
rz(-0.24793454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.068950653) q[0];
sx q[0];
rz(-2.729029) q[0];
sx q[0];
rz(2.5010342) q[0];
rz(-1.6259646) q[1];
sx q[1];
rz(-2.0104355) q[1];
sx q[1];
rz(-0.21387771) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2757077) q[0];
sx q[0];
rz(-1.2684221) q[0];
sx q[0];
rz(2.9024966) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1919695) q[2];
sx q[2];
rz(-2.2791499) q[2];
sx q[2];
rz(1.5713991) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7678309) q[1];
sx q[1];
rz(-1.6400178) q[1];
sx q[1];
rz(-0.94061416) q[1];
x q[2];
rz(1.6951896) q[3];
sx q[3];
rz(-0.52342192) q[3];
sx q[3];
rz(-1.0833797) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3512257) q[2];
sx q[2];
rz(-2.6595317) q[2];
sx q[2];
rz(-0.17769979) q[2];
rz(-1.7443582) q[3];
sx q[3];
rz(-1.9652365) q[3];
sx q[3];
rz(2.7241404) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.099667065) q[0];
sx q[0];
rz(-2.5302027) q[0];
sx q[0];
rz(2.0360816) q[0];
rz(0.28327495) q[1];
sx q[1];
rz(-2.6073644) q[1];
sx q[1];
rz(1.4615321) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.34237) q[0];
sx q[0];
rz(-0.91081753) q[0];
sx q[0];
rz(-0.8578542) q[0];
rz(2.2093532) q[2];
sx q[2];
rz(-1.9912212) q[2];
sx q[2];
rz(2.1115542) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.9859559) q[1];
sx q[1];
rz(-1.0105304) q[1];
sx q[1];
rz(1.7079279) q[1];
rz(-0.47847943) q[3];
sx q[3];
rz(-0.63204403) q[3];
sx q[3];
rz(1.5173591) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0854411) q[2];
sx q[2];
rz(-1.5458919) q[2];
sx q[2];
rz(2.936787) q[2];
rz(-1.0497931) q[3];
sx q[3];
rz(-2.864341) q[3];
sx q[3];
rz(-2.7753196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7751223) q[0];
sx q[0];
rz(-0.46638745) q[0];
sx q[0];
rz(-0.033576641) q[0];
rz(-1.4749984) q[1];
sx q[1];
rz(-0.74445236) q[1];
sx q[1];
rz(1.144145) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1468723) q[0];
sx q[0];
rz(-1.084274) q[0];
sx q[0];
rz(-1.5994607) q[0];
rz(-pi) q[1];
rz(1.8587684) q[2];
sx q[2];
rz(-2.8045553) q[2];
sx q[2];
rz(-2.2971643) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.7209789) q[1];
sx q[1];
rz(-0.79809626) q[1];
sx q[1];
rz(-2.6047203) q[1];
rz(-pi) q[2];
rz(-0.89544483) q[3];
sx q[3];
rz(-0.43206462) q[3];
sx q[3];
rz(2.2275782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.3860151) q[2];
sx q[2];
rz(-2.3365946) q[2];
sx q[2];
rz(1.2516652) q[2];
rz(1.3517316) q[3];
sx q[3];
rz(-0.91910619) q[3];
sx q[3];
rz(0.98021093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1397454) q[0];
sx q[0];
rz(-0.63848764) q[0];
sx q[0];
rz(1.8852604) q[0];
rz(-0.095257692) q[1];
sx q[1];
rz(-2.2525411) q[1];
sx q[1];
rz(2.6108066) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9448816) q[0];
sx q[0];
rz(-1.3232348) q[0];
sx q[0];
rz(2.976368) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.5321391) q[2];
sx q[2];
rz(-1.2529904) q[2];
sx q[2];
rz(-0.048887756) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3717613) q[1];
sx q[1];
rz(-2.8690845) q[1];
sx q[1];
rz(0.49401562) q[1];
rz(-0.23654273) q[3];
sx q[3];
rz(-0.45748392) q[3];
sx q[3];
rz(-1.2381697) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.329616) q[2];
sx q[2];
rz(-1.5831542) q[2];
sx q[2];
rz(2.6007593) q[2];
rz(-0.81816188) q[3];
sx q[3];
rz(-1.6988138) q[3];
sx q[3];
rz(2.5087859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5166017) q[0];
sx q[0];
rz(-1.7778968) q[0];
sx q[0];
rz(-2.7428108) q[0];
rz(-1.7575691) q[1];
sx q[1];
rz(-0.75474352) q[1];
sx q[1];
rz(3.0922281) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.931995) q[0];
sx q[0];
rz(-1.6315797) q[0];
sx q[0];
rz(1.6796735) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3081506) q[2];
sx q[2];
rz(-1.6742953) q[2];
sx q[2];
rz(-1.4348794) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.80222622) q[1];
sx q[1];
rz(-0.58103937) q[1];
sx q[1];
rz(1.4895579) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8667198) q[3];
sx q[3];
rz(-2.4913553) q[3];
sx q[3];
rz(0.45675983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.0177239) q[2];
sx q[2];
rz(-1.8327291) q[2];
sx q[2];
rz(-1.5105985) q[2];
rz(-0.92783582) q[3];
sx q[3];
rz(-1.8001013) q[3];
sx q[3];
rz(-1.9063037) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.31502003) q[0];
sx q[0];
rz(-2.6980504) q[0];
sx q[0];
rz(2.0499688) q[0];
rz(-1.0646959) q[1];
sx q[1];
rz(-0.5263435) q[1];
sx q[1];
rz(1.975504) q[1];
rz(0.76416107) q[2];
sx q[2];
rz(-0.66351009) q[2];
sx q[2];
rz(-1.4828975) q[2];
rz(2.8228685) q[3];
sx q[3];
rz(-2.095185) q[3];
sx q[3];
rz(-0.17879055) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
