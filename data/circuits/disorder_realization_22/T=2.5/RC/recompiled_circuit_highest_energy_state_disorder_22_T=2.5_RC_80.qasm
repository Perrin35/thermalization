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
rz(0.15777388) q[0];
sx q[0];
rz(-1.1717492) q[0];
sx q[0];
rz(1.2494614) q[0];
rz(-3.2818031) q[1];
sx q[1];
rz(1.6970716) q[1];
sx q[1];
rz(12.628218) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2551198) q[0];
sx q[0];
rz(-2.163871) q[0];
sx q[0];
rz(-2.5307239) q[0];
x q[1];
rz(-0.090422618) q[2];
sx q[2];
rz(-1.5711464) q[2];
sx q[2];
rz(1.4594913) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0272408) q[1];
sx q[1];
rz(-1.6682965) q[1];
sx q[1];
rz(-2.2932294) q[1];
rz(2.0560045) q[3];
sx q[3];
rz(-2.3923529) q[3];
sx q[3];
rz(-3.0390714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7141815) q[2];
sx q[2];
rz(-2.2990871) q[2];
sx q[2];
rz(0.25644914) q[2];
rz(-2.9668258) q[3];
sx q[3];
rz(-1.1656961) q[3];
sx q[3];
rz(-0.83785653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3278219) q[0];
sx q[0];
rz(-0.34532663) q[0];
sx q[0];
rz(-0.12369618) q[0];
rz(-2.9028614) q[1];
sx q[1];
rz(-0.78014603) q[1];
sx q[1];
rz(-3.0366268) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56120848) q[0];
sx q[0];
rz(-0.7155025) q[0];
sx q[0];
rz(-1.8010048) q[0];
x q[1];
rz(-2.6506422) q[2];
sx q[2];
rz(-0.46451312) q[2];
sx q[2];
rz(0.91769842) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.87890307) q[1];
sx q[1];
rz(-2.5715552) q[1];
sx q[1];
rz(-1.0703342) q[1];
rz(0.70521783) q[3];
sx q[3];
rz(-1.2949484) q[3];
sx q[3];
rz(0.57491702) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.3671942) q[2];
sx q[2];
rz(-1.4906733) q[2];
sx q[2];
rz(1.6801838) q[2];
rz(-1.2857619) q[3];
sx q[3];
rz(-1.6971842) q[3];
sx q[3];
rz(-0.27209601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9937781) q[0];
sx q[0];
rz(-0.21037924) q[0];
sx q[0];
rz(2.126597) q[0];
rz(0.89730942) q[1];
sx q[1];
rz(-1.4941447) q[1];
sx q[1];
rz(-0.6764594) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7583026) q[0];
sx q[0];
rz(-1.5967909) q[0];
sx q[0];
rz(2.1878408) q[0];
rz(-pi) q[1];
rz(-2.6482267) q[2];
sx q[2];
rz(-2.1378433) q[2];
sx q[2];
rz(0.21304785) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.8415422) q[1];
sx q[1];
rz(-1.9498697) q[1];
sx q[1];
rz(0.21684034) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5514938) q[3];
sx q[3];
rz(-2.3870644) q[3];
sx q[3];
rz(-1.7280662) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2066388) q[2];
sx q[2];
rz(-2.166344) q[2];
sx q[2];
rz(0.76033956) q[2];
rz(-1.2065411) q[3];
sx q[3];
rz(-0.81086603) q[3];
sx q[3];
rz(2.2981203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43561414) q[0];
sx q[0];
rz(-0.62788457) q[0];
sx q[0];
rz(1.5465558) q[0];
rz(-0.71890038) q[1];
sx q[1];
rz(-1.3201821) q[1];
sx q[1];
rz(-2.9100606) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0615094) q[0];
sx q[0];
rz(-0.69893796) q[0];
sx q[0];
rz(-0.7028347) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1377085) q[2];
sx q[2];
rz(-0.79266119) q[2];
sx q[2];
rz(-2.5156227) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(3.0445679) q[1];
sx q[1];
rz(-0.49355727) q[1];
sx q[1];
rz(2.943671) q[1];
x q[2];
rz(1.1778963) q[3];
sx q[3];
rz(-2.1390954) q[3];
sx q[3];
rz(1.5400122) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.71451688) q[2];
sx q[2];
rz(-2.7526553) q[2];
sx q[2];
rz(2.5330949) q[2];
rz(0.19615873) q[3];
sx q[3];
rz(-1.4213296) q[3];
sx q[3];
rz(2.1430446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(1.8683559) q[0];
sx q[0];
rz(-0.65160692) q[0];
sx q[0];
rz(2.3625145) q[0];
rz(-0.92998663) q[1];
sx q[1];
rz(-1.9735186) q[1];
sx q[1];
rz(-1.9680061) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1176276) q[0];
sx q[0];
rz(-1.574734) q[0];
sx q[0];
rz(-1.2327475) q[0];
rz(-pi) q[1];
x q[1];
rz(0.13466882) q[2];
sx q[2];
rz(-2.1818518) q[2];
sx q[2];
rz(-2.5650283) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.47347676) q[1];
sx q[1];
rz(-1.887823) q[1];
sx q[1];
rz(-3.0739529) q[1];
rz(0.39802528) q[3];
sx q[3];
rz(-2.6809691) q[3];
sx q[3];
rz(-1.6470136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.32213) q[2];
sx q[2];
rz(-2.9746015) q[2];
sx q[2];
rz(0.67506153) q[2];
rz(-2.2760462) q[3];
sx q[3];
rz(-1.5138488) q[3];
sx q[3];
rz(-1.9338098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.753767) q[0];
sx q[0];
rz(-0.017711552) q[0];
sx q[0];
rz(-2.2702763) q[0];
rz(0.21698347) q[1];
sx q[1];
rz(-1.5419518) q[1];
sx q[1];
rz(-1.8035696) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8080475) q[0];
sx q[0];
rz(-1.4879972) q[0];
sx q[0];
rz(-2.4890578) q[0];
rz(-pi) q[1];
rz(1.4124521) q[2];
sx q[2];
rz(-0.73190231) q[2];
sx q[2];
rz(-0.5754234) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8585526) q[1];
sx q[1];
rz(-1.5225019) q[1];
sx q[1];
rz(-1.7142536) q[1];
x q[2];
rz(1.7609672) q[3];
sx q[3];
rz(-2.0479879) q[3];
sx q[3];
rz(2.793345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-3.1285105) q[2];
sx q[2];
rz(-1.4980114) q[2];
sx q[2];
rz(-0.80005542) q[2];
rz(-0.99177805) q[3];
sx q[3];
rz(-2.1938775) q[3];
sx q[3];
rz(2.1984524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-0.25471383) q[0];
sx q[0];
rz(-0.44097057) q[0];
sx q[0];
rz(0.15175858) q[0];
rz(1.4166547) q[1];
sx q[1];
rz(-2.2817426) q[1];
sx q[1];
rz(-0.83622611) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4164616) q[0];
sx q[0];
rz(-2.2380479) q[0];
sx q[0];
rz(0.50531549) q[0];
x q[1];
rz(2.5649037) q[2];
sx q[2];
rz(-1.5637914) q[2];
sx q[2];
rz(3.0617065) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6265833) q[1];
sx q[1];
rz(-1.3224416) q[1];
sx q[1];
rz(-0.088661389) q[1];
rz(-2.6760929) q[3];
sx q[3];
rz(-1.4537849) q[3];
sx q[3];
rz(2.4586364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.4862711) q[2];
sx q[2];
rz(-0.063491193) q[2];
sx q[2];
rz(3.0625694) q[2];
rz(-2.4231353) q[3];
sx q[3];
rz(-1.5114762) q[3];
sx q[3];
rz(-0.92459905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.786161) q[0];
sx q[0];
rz(-2.5434255) q[0];
sx q[0];
rz(0.86791903) q[0];
rz(2.4527841) q[1];
sx q[1];
rz(-2.8354366) q[1];
sx q[1];
rz(0.93592962) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.82831885) q[0];
sx q[0];
rz(-2.5848645) q[0];
sx q[0];
rz(0.3484881) q[0];
rz(0.88603381) q[2];
sx q[2];
rz(-2.6305969) q[2];
sx q[2];
rz(0.053879468) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.28081255) q[1];
sx q[1];
rz(-1.6242923) q[1];
sx q[1];
rz(-1.2455275) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5722365) q[3];
sx q[3];
rz(-2.3142696) q[3];
sx q[3];
rz(-1.8515406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1615289) q[2];
sx q[2];
rz(-2.4825725) q[2];
sx q[2];
rz(-2.8847983) q[2];
rz(-1.0181381) q[3];
sx q[3];
rz(-2.0223821) q[3];
sx q[3];
rz(0.97091466) q[3];
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
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37860206) q[0];
sx q[0];
rz(-2.584223) q[0];
sx q[0];
rz(-0.85365224) q[0];
rz(-1.8866106) q[1];
sx q[1];
rz(-1.1458784) q[1];
sx q[1];
rz(0.34057158) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.69215323) q[0];
sx q[0];
rz(-1.7652349) q[0];
sx q[0];
rz(2.9783335) q[0];
x q[1];
rz(-2.8355424) q[2];
sx q[2];
rz(-2.4444164) q[2];
sx q[2];
rz(-2.9070284) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7813959) q[1];
sx q[1];
rz(-0.92222795) q[1];
sx q[1];
rz(0.99332033) q[1];
rz(-pi) q[2];
rz(2.0331618) q[3];
sx q[3];
rz(-1.3232854) q[3];
sx q[3];
rz(-0.68456291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.42369947) q[2];
sx q[2];
rz(-1.3807978) q[2];
sx q[2];
rz(-2.763486) q[2];
rz(0.2374436) q[3];
sx q[3];
rz(-2.0707371) q[3];
sx q[3];
rz(2.8606991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0963999) q[0];
sx q[0];
rz(-2.0572331) q[0];
sx q[0];
rz(0.64724809) q[0];
rz(-1.1680394) q[1];
sx q[1];
rz(-0.85473514) q[1];
sx q[1];
rz(-2.7580269) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3824359) q[0];
sx q[0];
rz(-1.54017) q[0];
sx q[0];
rz(-1.5923862) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6099989) q[2];
sx q[2];
rz(-1.2832912) q[2];
sx q[2];
rz(2.2026041) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.3473222) q[1];
sx q[1];
rz(-1.3065265) q[1];
sx q[1];
rz(1.3572973) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0342919) q[3];
sx q[3];
rz(-1.0578053) q[3];
sx q[3];
rz(1.9657621) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8770404) q[2];
sx q[2];
rz(-0.88901797) q[2];
sx q[2];
rz(1.5569347) q[2];
rz(1.5133739) q[3];
sx q[3];
rz(-1.7729365) q[3];
sx q[3];
rz(0.42678601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3601111) q[0];
sx q[0];
rz(-1.7088912) q[0];
sx q[0];
rz(0.29722469) q[0];
rz(-1.9602641) q[1];
sx q[1];
rz(-1.344463) q[1];
sx q[1];
rz(-0.32483473) q[1];
rz(1.0931821) q[2];
sx q[2];
rz(-1.3543159) q[2];
sx q[2];
rz(-1.3456624) q[2];
rz(-1.539113) q[3];
sx q[3];
rz(-2.0551695) q[3];
sx q[3];
rz(2.7885128) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
