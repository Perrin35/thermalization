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
rz(0.76437104) q[0];
sx q[0];
rz(-1.3410913) q[0];
sx q[0];
rz(2.2361225) q[0];
rz(1.5068997) q[1];
sx q[1];
rz(-2.1807179) q[1];
sx q[1];
rz(-1.6406055) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1849821) q[0];
sx q[0];
rz(-2.2652103) q[0];
sx q[0];
rz(2.882769) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9328961) q[2];
sx q[2];
rz(-2.675867) q[2];
sx q[2];
rz(3.1052239) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.49351827) q[1];
sx q[1];
rz(-1.5941125) q[1];
sx q[1];
rz(0.042453153) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1330058) q[3];
sx q[3];
rz(-1.2744181) q[3];
sx q[3];
rz(-2.5135771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.91244873) q[2];
sx q[2];
rz(-1.7488166) q[2];
sx q[2];
rz(-0.32192117) q[2];
rz(-2.1866482) q[3];
sx q[3];
rz(-0.82986444) q[3];
sx q[3];
rz(-1.1436852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2442653) q[0];
sx q[0];
rz(-1.0538415) q[0];
sx q[0];
rz(0.51189297) q[0];
rz(1.5882675) q[1];
sx q[1];
rz(-0.48738185) q[1];
sx q[1];
rz(1.1221251) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7417541) q[0];
sx q[0];
rz(-1.4692093) q[0];
sx q[0];
rz(-0.31272696) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1013158) q[2];
sx q[2];
rz(-1.1415408) q[2];
sx q[2];
rz(-1.5778936) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.304804) q[1];
sx q[1];
rz(-1.2366017) q[1];
sx q[1];
rz(2.0242974) q[1];
rz(-2.4347794) q[3];
sx q[3];
rz(-1.9403096) q[3];
sx q[3];
rz(0.57716864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.7684218) q[2];
sx q[2];
rz(-1.3139407) q[2];
sx q[2];
rz(2.6785417) q[2];
rz(-0.57524663) q[3];
sx q[3];
rz(-1.3653711) q[3];
sx q[3];
rz(-3.0423394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7333882) q[0];
sx q[0];
rz(-0.84771228) q[0];
sx q[0];
rz(2.7114765) q[0];
rz(0.46562132) q[1];
sx q[1];
rz(-0.72048134) q[1];
sx q[1];
rz(-2.5045085) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.67733902) q[0];
sx q[0];
rz(-1.5298944) q[0];
sx q[0];
rz(-3.1275355) q[0];
rz(1.5970206) q[2];
sx q[2];
rz(-2.3560212) q[2];
sx q[2];
rz(0.77411133) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.068598572) q[1];
sx q[1];
rz(-1.322876) q[1];
sx q[1];
rz(-0.74241717) q[1];
rz(-1.1807084) q[3];
sx q[3];
rz(-1.0849139) q[3];
sx q[3];
rz(2.7220243) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.78274) q[2];
sx q[2];
rz(-1.0461067) q[2];
sx q[2];
rz(2.4165238) q[2];
rz(1.8105761) q[3];
sx q[3];
rz(-2.1264117) q[3];
sx q[3];
rz(-0.63841188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8722039) q[0];
sx q[0];
rz(-1.6870455) q[0];
sx q[0];
rz(1.697502) q[0];
rz(-3.0929502) q[1];
sx q[1];
rz(-1.8154058) q[1];
sx q[1];
rz(2.8526502) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64000696) q[0];
sx q[0];
rz(-1.3292392) q[0];
sx q[0];
rz(1.1904741) q[0];
rz(1.7554531) q[2];
sx q[2];
rz(-0.15310213) q[2];
sx q[2];
rz(-2.7775922) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.1565981) q[1];
sx q[1];
rz(-0.20594507) q[1];
sx q[1];
rz(2.1316281) q[1];
rz(-1.5163317) q[3];
sx q[3];
rz(-0.90622444) q[3];
sx q[3];
rz(1.8502473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.85663969) q[2];
sx q[2];
rz(-0.53264561) q[2];
sx q[2];
rz(0.3802158) q[2];
rz(0.63816655) q[3];
sx q[3];
rz(-1.2297945) q[3];
sx q[3];
rz(-1.1285454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83288348) q[0];
sx q[0];
rz(-0.55279151) q[0];
sx q[0];
rz(1.0943476) q[0];
rz(0.87567466) q[1];
sx q[1];
rz(-0.62962571) q[1];
sx q[1];
rz(2.0424776) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7849784) q[0];
sx q[0];
rz(-1.4334502) q[0];
sx q[0];
rz(-0.38909586) q[0];
rz(-pi) q[1];
rz(2.3540456) q[2];
sx q[2];
rz(-0.468245) q[2];
sx q[2];
rz(-2.5494273) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-3.0620516) q[1];
sx q[1];
rz(-1.4126443) q[1];
sx q[1];
rz(-1.016721) q[1];
rz(-pi) q[2];
x q[2];
rz(0.15954475) q[3];
sx q[3];
rz(-2.2189848) q[3];
sx q[3];
rz(-1.2153347) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2431474) q[2];
sx q[2];
rz(-1.7869608) q[2];
sx q[2];
rz(-2.1649427) q[2];
rz(0.53168932) q[3];
sx q[3];
rz(-0.44638005) q[3];
sx q[3];
rz(-0.71438742) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0831182) q[0];
sx q[0];
rz(-1.1019022) q[0];
sx q[0];
rz(1.8699159) q[0];
rz(2.7187128) q[1];
sx q[1];
rz(-1.5300749) q[1];
sx q[1];
rz(1.6082825) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6278224) q[0];
sx q[0];
rz(-1.5585527) q[0];
sx q[0];
rz(-1.5041385) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6688231) q[2];
sx q[2];
rz(-1.5614126) q[2];
sx q[2];
rz(-0.90017747) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.21580885) q[1];
sx q[1];
rz(-1.5152351) q[1];
sx q[1];
rz(2.8381398) q[1];
rz(-pi) q[2];
x q[2];
rz(0.23171429) q[3];
sx q[3];
rz(-0.82666708) q[3];
sx q[3];
rz(0.8207013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.61234683) q[2];
sx q[2];
rz(-1.2924478) q[2];
sx q[2];
rz(-2.7563654) q[2];
rz(0.084608229) q[3];
sx q[3];
rz(-2.707983) q[3];
sx q[3];
rz(2.2658074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.92152921) q[0];
sx q[0];
rz(-1.055287) q[0];
sx q[0];
rz(-0.39749843) q[0];
rz(2.6037604) q[1];
sx q[1];
rz(-0.42306867) q[1];
sx q[1];
rz(-0.0040815512) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57848141) q[0];
sx q[0];
rz(-2.3627776) q[0];
sx q[0];
rz(-1.388873) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.053431) q[2];
sx q[2];
rz(-0.74649278) q[2];
sx q[2];
rz(1.5160402) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.38975484) q[1];
sx q[1];
rz(-2.1389277) q[1];
sx q[1];
rz(2.0271432) q[1];
x q[2];
rz(2.2479731) q[3];
sx q[3];
rz(-2.1359518) q[3];
sx q[3];
rz(-2.2806185) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.5092545) q[2];
sx q[2];
rz(-1.961144) q[2];
sx q[2];
rz(-2.9285367) q[2];
rz(-2.3325855) q[3];
sx q[3];
rz(-0.041497858) q[3];
sx q[3];
rz(-1.8607148) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0179366) q[0];
sx q[0];
rz(-1.2024711) q[0];
sx q[0];
rz(1.9785471) q[0];
rz(-0.2991547) q[1];
sx q[1];
rz(-1.2048293) q[1];
sx q[1];
rz(0.97602731) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.19071391) q[0];
sx q[0];
rz(-1.4602721) q[0];
sx q[0];
rz(-2.8948363) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.15248044) q[2];
sx q[2];
rz(-1.7465542) q[2];
sx q[2];
rz(-0.7712785) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2175156) q[1];
sx q[1];
rz(-1.2189606) q[1];
sx q[1];
rz(2.887243) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7175432) q[3];
sx q[3];
rz(-1.3433787) q[3];
sx q[3];
rz(-3.1229179) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.8396847) q[2];
sx q[2];
rz(-2.8161616) q[2];
sx q[2];
rz(3.0111266) q[2];
rz(-1.3236375) q[3];
sx q[3];
rz(-1.2490844) q[3];
sx q[3];
rz(-0.29449335) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94006222) q[0];
sx q[0];
rz(-0.37684965) q[0];
sx q[0];
rz(-1.4991722) q[0];
rz(2.6014853) q[1];
sx q[1];
rz(-0.79743782) q[1];
sx q[1];
rz(1.3444208) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.62468707) q[0];
sx q[0];
rz(-0.4962099) q[0];
sx q[0];
rz(2.2115256) q[0];
rz(-pi) q[1];
x q[1];
rz(0.62144582) q[2];
sx q[2];
rz(-1.1100195) q[2];
sx q[2];
rz(-0.64475343) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0186179) q[1];
sx q[1];
rz(-2.5135165) q[1];
sx q[1];
rz(1.1372034) q[1];
rz(-pi) q[2];
rz(-2.5193941) q[3];
sx q[3];
rz(-0.95131627) q[3];
sx q[3];
rz(3.0752237) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.53238791) q[2];
sx q[2];
rz(-2.8049073) q[2];
sx q[2];
rz(1.6877635) q[2];
rz(-2.7135571) q[3];
sx q[3];
rz(-1.5902218) q[3];
sx q[3];
rz(1.154703) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5880599) q[0];
sx q[0];
rz(-0.45192161) q[0];
sx q[0];
rz(1.4916627) q[0];
rz(-0.55117575) q[1];
sx q[1];
rz(-2.1291514) q[1];
sx q[1];
rz(2.5490882) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6766178) q[0];
sx q[0];
rz(-1.3335557) q[0];
sx q[0];
rz(0.28740164) q[0];
rz(-0.5430605) q[2];
sx q[2];
rz(-2.0351699) q[2];
sx q[2];
rz(-0.65242243) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.77280918) q[1];
sx q[1];
rz(-1.9516489) q[1];
sx q[1];
rz(1.8541149) q[1];
x q[2];
rz(-3.0499712) q[3];
sx q[3];
rz(-2.6406248) q[3];
sx q[3];
rz(-0.17551455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.7117915) q[2];
sx q[2];
rz(-1.0415123) q[2];
sx q[2];
rz(0.05833021) q[2];
rz(-0.37060261) q[3];
sx q[3];
rz(-2.8648418) q[3];
sx q[3];
rz(-3.1378194) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2881099) q[0];
sx q[0];
rz(-1.7848889) q[0];
sx q[0];
rz(-3.0369192) q[0];
rz(-1.3660322) q[1];
sx q[1];
rz(-0.80148253) q[1];
sx q[1];
rz(-2.186224) q[1];
rz(-1.3132533) q[2];
sx q[2];
rz(-1.0124442) q[2];
sx q[2];
rz(-0.53895216) q[2];
rz(-0.80647918) q[3];
sx q[3];
rz(-2.4462593) q[3];
sx q[3];
rz(-1.367955) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
