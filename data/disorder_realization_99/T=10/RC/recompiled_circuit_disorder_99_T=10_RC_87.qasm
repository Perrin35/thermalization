OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.0302304) q[0];
sx q[0];
rz(-0.65519873) q[0];
sx q[0];
rz(0.97638786) q[0];
rz(-0.916565) q[1];
sx q[1];
rz(-0.7601127) q[1];
sx q[1];
rz(2.3488933) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22098955) q[0];
sx q[0];
rz(-1.594829) q[0];
sx q[0];
rz(-2.241797) q[0];
x q[1];
rz(2.5976074) q[2];
sx q[2];
rz(-1.7684801) q[2];
sx q[2];
rz(-2.939784) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.77289171) q[1];
sx q[1];
rz(-1.643631) q[1];
sx q[1];
rz(-1.1198977) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.84150746) q[3];
sx q[3];
rz(-0.91472799) q[3];
sx q[3];
rz(-0.63667242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.67316002) q[2];
sx q[2];
rz(-2.1940553) q[2];
sx q[2];
rz(0.49896487) q[2];
rz(-2.7178102) q[3];
sx q[3];
rz(-1.2200004) q[3];
sx q[3];
rz(3.1288778) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38133165) q[0];
sx q[0];
rz(-1.6853764) q[0];
sx q[0];
rz(-0.99200845) q[0];
rz(-0.17653067) q[1];
sx q[1];
rz(-2.9361528) q[1];
sx q[1];
rz(-2.7761249) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.394709) q[0];
sx q[0];
rz(-2.611126) q[0];
sx q[0];
rz(-0.51325004) q[0];
rz(-pi) q[1];
rz(0.59132691) q[2];
sx q[2];
rz(-2.1543243) q[2];
sx q[2];
rz(3.025625) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6668015) q[1];
sx q[1];
rz(-2.431369) q[1];
sx q[1];
rz(0.53479654) q[1];
rz(-pi) q[2];
x q[2];
rz(2.9404853) q[3];
sx q[3];
rz(-1.6893941) q[3];
sx q[3];
rz(-3.0191641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.83313292) q[2];
sx q[2];
rz(-0.59811991) q[2];
sx q[2];
rz(-1.0773405) q[2];
rz(-0.16263738) q[3];
sx q[3];
rz(-1.2626303) q[3];
sx q[3];
rz(-1.8796273) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4804374) q[0];
sx q[0];
rz(-0.85395956) q[0];
sx q[0];
rz(-0.61082947) q[0];
rz(-1.707533) q[1];
sx q[1];
rz(-1.779498) q[1];
sx q[1];
rz(2.6909713) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.075899344) q[0];
sx q[0];
rz(-1.6487299) q[0];
sx q[0];
rz(-3.1083641) q[0];
rz(-pi) q[1];
rz(-2.8350416) q[2];
sx q[2];
rz(-1.2129285) q[2];
sx q[2];
rz(-1.7265665) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.9420942) q[1];
sx q[1];
rz(-0.71951413) q[1];
sx q[1];
rz(2.662563) q[1];
rz(-pi) q[2];
rz(-1.2080473) q[3];
sx q[3];
rz(-2.4137073) q[3];
sx q[3];
rz(-2.4129652) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0314363) q[2];
sx q[2];
rz(-2.0560196) q[2];
sx q[2];
rz(-1.191656) q[2];
rz(-2.2971161) q[3];
sx q[3];
rz(-2.8561487) q[3];
sx q[3];
rz(0.58117956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.095117005) q[0];
sx q[0];
rz(-1.3866383) q[0];
sx q[0];
rz(2.0003831) q[0];
rz(1.8449239) q[1];
sx q[1];
rz(-0.43162391) q[1];
sx q[1];
rz(0.30532125) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5019708) q[0];
sx q[0];
rz(-2.1272215) q[0];
sx q[0];
rz(-1.9407546) q[0];
rz(-pi) q[1];
x q[1];
rz(2.4112169) q[2];
sx q[2];
rz(-2.5502898) q[2];
sx q[2];
rz(3.0808466) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.2480337) q[1];
sx q[1];
rz(-2.3602924) q[1];
sx q[1];
rz(1.7209189) q[1];
x q[2];
rz(-2.3577945) q[3];
sx q[3];
rz(-0.67516203) q[3];
sx q[3];
rz(2.6019707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.9332463) q[2];
sx q[2];
rz(-1.5657921) q[2];
sx q[2];
rz(-0.34710458) q[2];
rz(0.38337213) q[3];
sx q[3];
rz(-2.255286) q[3];
sx q[3];
rz(-2.0736407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40889302) q[0];
sx q[0];
rz(-2.2786031) q[0];
sx q[0];
rz(2.6026759) q[0];
rz(-1.7319038) q[1];
sx q[1];
rz(-0.46202818) q[1];
sx q[1];
rz(2.7358823) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.51605326) q[0];
sx q[0];
rz(-2.2279983) q[0];
sx q[0];
rz(-1.9198734) q[0];
rz(-pi) q[1];
rz(-1.97498) q[2];
sx q[2];
rz(-2.5378072) q[2];
sx q[2];
rz(-3.0886138) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.87863641) q[1];
sx q[1];
rz(-2.7107781) q[1];
sx q[1];
rz(-0.93114914) q[1];
x q[2];
rz(1.3375072) q[3];
sx q[3];
rz(-0.81987112) q[3];
sx q[3];
rz(-1.6398721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.47200176) q[2];
sx q[2];
rz(-0.80299032) q[2];
sx q[2];
rz(1.6873138) q[2];
rz(-0.83868319) q[3];
sx q[3];
rz(-1.8554153) q[3];
sx q[3];
rz(1.3735166) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8307761) q[0];
sx q[0];
rz(-2.5287479) q[0];
sx q[0];
rz(2.7344761) q[0];
rz(-2.282417) q[1];
sx q[1];
rz(-1.5770864) q[1];
sx q[1];
rz(-1.3173332) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.28535298) q[0];
sx q[0];
rz(-2.1137153) q[0];
sx q[0];
rz(3.1154484) q[0];
x q[1];
rz(0.66770422) q[2];
sx q[2];
rz(-1.0945079) q[2];
sx q[2];
rz(-2.8712733) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.6588905) q[1];
sx q[1];
rz(-1.5679949) q[1];
sx q[1];
rz(-0.83562327) q[1];
x q[2];
rz(2.7139211) q[3];
sx q[3];
rz(-2.1139522) q[3];
sx q[3];
rz(-3.1205408) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.09981) q[2];
sx q[2];
rz(-0.87974292) q[2];
sx q[2];
rz(0.89522925) q[2];
rz(1.9334531) q[3];
sx q[3];
rz(-0.67966214) q[3];
sx q[3];
rz(-2.4334548) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5642501) q[0];
sx q[0];
rz(-1.0471434) q[0];
sx q[0];
rz(2.0136925) q[0];
rz(-0.212542) q[1];
sx q[1];
rz(-1.3393341) q[1];
sx q[1];
rz(1.2316661) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2736149) q[0];
sx q[0];
rz(-0.6379188) q[0];
sx q[0];
rz(0.62423737) q[0];
rz(-2.3627794) q[2];
sx q[2];
rz(-0.80873571) q[2];
sx q[2];
rz(0.29701172) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.61464804) q[1];
sx q[1];
rz(-1.467418) q[1];
sx q[1];
rz(2.9109216) q[1];
x q[2];
rz(1.6024186) q[3];
sx q[3];
rz(-0.60472371) q[3];
sx q[3];
rz(-1.5865758) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3125375) q[2];
sx q[2];
rz(-2.2237491) q[2];
sx q[2];
rz(0.60116872) q[2];
rz(0.51856315) q[3];
sx q[3];
rz(-2.8278567) q[3];
sx q[3];
rz(1.7140088) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1475875) q[0];
sx q[0];
rz(-3.1412791) q[0];
sx q[0];
rz(3.0684379) q[0];
rz(2.1144497) q[1];
sx q[1];
rz(-1.606753) q[1];
sx q[1];
rz(-2.5491319) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.52623977) q[0];
sx q[0];
rz(-2.549299) q[0];
sx q[0];
rz(-0.083239716) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4249435) q[2];
sx q[2];
rz(-0.71084329) q[2];
sx q[2];
rz(-2.5568642) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2295099) q[1];
sx q[1];
rz(-0.76369897) q[1];
sx q[1];
rz(-2.0565226) q[1];
rz(2.8665364) q[3];
sx q[3];
rz(-1.1797136) q[3];
sx q[3];
rz(-2.070379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.5804194) q[2];
sx q[2];
rz(-1.6757123) q[2];
sx q[2];
rz(-0.44211659) q[2];
rz(0.90028611) q[3];
sx q[3];
rz(-1.0628275) q[3];
sx q[3];
rz(0.013307868) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1135547) q[0];
sx q[0];
rz(-2.2224764) q[0];
sx q[0];
rz(-1.8369209) q[0];
rz(-1.4470944) q[1];
sx q[1];
rz(-1.1728975) q[1];
sx q[1];
rz(-0.56217271) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5484555) q[0];
sx q[0];
rz(-2.2459185) q[0];
sx q[0];
rz(1.1086585) q[0];
rz(-pi) q[1];
rz(1.4116889) q[2];
sx q[2];
rz(-0.53630398) q[2];
sx q[2];
rz(-1.916534) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.11645424) q[1];
sx q[1];
rz(-2.0471731) q[1];
sx q[1];
rz(-0.78671793) q[1];
rz(-pi) q[2];
rz(-2.3638944) q[3];
sx q[3];
rz(-0.54293984) q[3];
sx q[3];
rz(2.267595) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8582981) q[2];
sx q[2];
rz(-1.4026182) q[2];
sx q[2];
rz(0.98158681) q[2];
rz(-3.1372519) q[3];
sx q[3];
rz(-2.0948295) q[3];
sx q[3];
rz(2.5481352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3146064) q[0];
sx q[0];
rz(-1.9208603) q[0];
sx q[0];
rz(-3.0944968) q[0];
rz(1.3866407) q[1];
sx q[1];
rz(-1.9853233) q[1];
sx q[1];
rz(-0.047853619) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3210147) q[0];
sx q[0];
rz(-0.65549675) q[0];
sx q[0];
rz(1.5232248) q[0];
rz(-pi) q[1];
rz(-1.1282975) q[2];
sx q[2];
rz(-1.6167165) q[2];
sx q[2];
rz(-2.0890582) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.3198941) q[1];
sx q[1];
rz(-2.387429) q[1];
sx q[1];
rz(2.2241705) q[1];
rz(-pi) q[2];
rz(3.0129274) q[3];
sx q[3];
rz(-0.4056969) q[3];
sx q[3];
rz(-2.2588244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4973267) q[2];
sx q[2];
rz(-1.9571807) q[2];
sx q[2];
rz(-2.2052374) q[2];
rz(-2.3406773) q[3];
sx q[3];
rz(-0.51910669) q[3];
sx q[3];
rz(0.67805964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3180278) q[0];
sx q[0];
rz(-0.93946539) q[0];
sx q[0];
rz(-2.9325624) q[0];
rz(-2.6387852) q[1];
sx q[1];
rz(-1.8223169) q[1];
sx q[1];
rz(1.0261789) q[1];
rz(0.55124333) q[2];
sx q[2];
rz(-1.5894645) q[2];
sx q[2];
rz(-0.20655256) q[2];
rz(-2.7754178) q[3];
sx q[3];
rz(-1.1911285) q[3];
sx q[3];
rz(1.7366684) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
