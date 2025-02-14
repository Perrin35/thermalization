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
rz(-0.40773243) q[0];
sx q[0];
rz(-1.9324349) q[0];
sx q[0];
rz(1.6785167) q[0];
rz(-2.8830124) q[1];
sx q[1];
rz(-1.2376031) q[1];
sx q[1];
rz(0.13374506) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8834849) q[0];
sx q[0];
rz(-2.9843669) q[0];
sx q[0];
rz(1.5327461) q[0];
rz(-pi) q[1];
x q[1];
rz(0.76579801) q[2];
sx q[2];
rz(-2.4054689) q[2];
sx q[2];
rz(-2.78139) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.3901538) q[1];
sx q[1];
rz(-0.62733106) q[1];
sx q[1];
rz(2.6120758) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.7236563) q[3];
sx q[3];
rz(-1.5709849) q[3];
sx q[3];
rz(2.3589695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3818843) q[2];
sx q[2];
rz(-2.3298161) q[2];
sx q[2];
rz(-1.3105357) q[2];
rz(-1.367502) q[3];
sx q[3];
rz(-2.01367) q[3];
sx q[3];
rz(3.1254356) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1644208) q[0];
sx q[0];
rz(-1.9855969) q[0];
sx q[0];
rz(-1.0154065) q[0];
rz(0.39335355) q[1];
sx q[1];
rz(-1.669408) q[1];
sx q[1];
rz(-0.70158395) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6664826) q[0];
sx q[0];
rz(-2.283264) q[0];
sx q[0];
rz(3.0300167) q[0];
rz(-pi) q[1];
rz(0.010741269) q[2];
sx q[2];
rz(-1.2678384) q[2];
sx q[2];
rz(-3.0634653) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.1893411) q[1];
sx q[1];
rz(-2.8165112) q[1];
sx q[1];
rz(1.3334683) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.75922482) q[3];
sx q[3];
rz(-2.1442778) q[3];
sx q[3];
rz(-0.92347565) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.027448805) q[2];
sx q[2];
rz(-1.0237209) q[2];
sx q[2];
rz(1.1055498) q[2];
rz(-0.28087273) q[3];
sx q[3];
rz(-0.61087817) q[3];
sx q[3];
rz(-0.15225473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1059145) q[0];
sx q[0];
rz(-2.4690101) q[0];
sx q[0];
rz(-1.0726844) q[0];
rz(-0.013280344) q[1];
sx q[1];
rz(-1.599879) q[1];
sx q[1];
rz(0.57659155) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.428513) q[0];
sx q[0];
rz(-1.858485) q[0];
sx q[0];
rz(1.6384533) q[0];
rz(-1.7037839) q[2];
sx q[2];
rz(-1.4879359) q[2];
sx q[2];
rz(1.8216009) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.62200786) q[1];
sx q[1];
rz(-1.7652006) q[1];
sx q[1];
rz(2.0561051) q[1];
x q[2];
rz(2.1284038) q[3];
sx q[3];
rz(-1.3609145) q[3];
sx q[3];
rz(-2.5887973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.85792929) q[2];
sx q[2];
rz(-0.89202213) q[2];
sx q[2];
rz(-0.64209765) q[2];
rz(0.56504956) q[3];
sx q[3];
rz(-2.8665906) q[3];
sx q[3];
rz(-0.16998418) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
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
rz(0.86244407) q[0];
sx q[0];
rz(-1.886241) q[0];
sx q[0];
rz(2.8992262) q[0];
rz(0.53388059) q[1];
sx q[1];
rz(-2.4568074) q[1];
sx q[1];
rz(1.4885611) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87391728) q[0];
sx q[0];
rz(-0.59989625) q[0];
sx q[0];
rz(-3.0156662) q[0];
rz(-pi) q[1];
x q[1];
rz(2.7127391) q[2];
sx q[2];
rz(-1.2786713) q[2];
sx q[2];
rz(-3.0512864) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.44806717) q[1];
sx q[1];
rz(-2.2972745) q[1];
sx q[1];
rz(-1.179715) q[1];
x q[2];
rz(2.8418301) q[3];
sx q[3];
rz(-1.5325755) q[3];
sx q[3];
rz(-0.37118566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.7786467) q[2];
sx q[2];
rz(-1.2641509) q[2];
sx q[2];
rz(-1.470559) q[2];
rz(-0.10739022) q[3];
sx q[3];
rz(-1.3067747) q[3];
sx q[3];
rz(-1.86514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.030468) q[0];
sx q[0];
rz(-2.680439) q[0];
sx q[0];
rz(1.2029458) q[0];
rz(2.2466834) q[1];
sx q[1];
rz(-1.1824965) q[1];
sx q[1];
rz(-0.21200171) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.934904) q[0];
sx q[0];
rz(-2.3634533) q[0];
sx q[0];
rz(-2.2147793) q[0];
rz(-pi) q[1];
x q[1];
rz(0.49290979) q[2];
sx q[2];
rz(-2.5625336) q[2];
sx q[2];
rz(-2.2185203) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.70658503) q[1];
sx q[1];
rz(-2.3870578) q[1];
sx q[1];
rz(-0.88891397) q[1];
rz(1.5552484) q[3];
sx q[3];
rz(-2.5882396) q[3];
sx q[3];
rz(-1.4912332) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2682858) q[2];
sx q[2];
rz(-0.77750677) q[2];
sx q[2];
rz(0.11534616) q[2];
rz(-2.1558732) q[3];
sx q[3];
rz(-1.7306354) q[3];
sx q[3];
rz(0.43578291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91144052) q[0];
sx q[0];
rz(-1.6807012) q[0];
sx q[0];
rz(0.69754115) q[0];
rz(-2.7279834) q[1];
sx q[1];
rz(-1.4613084) q[1];
sx q[1];
rz(-2.4381309) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7005806) q[0];
sx q[0];
rz(-1.3003775) q[0];
sx q[0];
rz(-1.7187814) q[0];
x q[1];
rz(-2.1036316) q[2];
sx q[2];
rz(-2.46358) q[2];
sx q[2];
rz(0.17282669) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.9338434) q[1];
sx q[1];
rz(-1.9393171) q[1];
sx q[1];
rz(-0.50362419) q[1];
x q[2];
rz(-1.2715152) q[3];
sx q[3];
rz(-1.3320775) q[3];
sx q[3];
rz(-1.6629926) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.10282639) q[2];
sx q[2];
rz(-0.26831728) q[2];
sx q[2];
rz(1.212567) q[2];
rz(0.28193998) q[3];
sx q[3];
rz(-1.5423256) q[3];
sx q[3];
rz(-2.6350002) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3459699) q[0];
sx q[0];
rz(-2.3891734) q[0];
sx q[0];
rz(0.84683007) q[0];
rz(-2.1961191) q[1];
sx q[1];
rz(-1.5676326) q[1];
sx q[1];
rz(-2.4623154) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4436809) q[0];
sx q[0];
rz(-2.6583932) q[0];
sx q[0];
rz(2.9761821) q[0];
x q[1];
rz(0.62394721) q[2];
sx q[2];
rz(-1.3610164) q[2];
sx q[2];
rz(-2.7854475) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5838562) q[1];
sx q[1];
rz(-2.5917555) q[1];
sx q[1];
rz(-1.932895) q[1];
rz(1.6369989) q[3];
sx q[3];
rz(-0.28488628) q[3];
sx q[3];
rz(-1.5004304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.2519553) q[2];
sx q[2];
rz(-2.1869982) q[2];
sx q[2];
rz(1.4944705) q[2];
rz(2.5877118) q[3];
sx q[3];
rz(-1.5364105) q[3];
sx q[3];
rz(1.7249829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81903356) q[0];
sx q[0];
rz(-0.77545866) q[0];
sx q[0];
rz(-1.092528) q[0];
rz(0.93387261) q[1];
sx q[1];
rz(-1.7364419) q[1];
sx q[1];
rz(0.78528231) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99150554) q[0];
sx q[0];
rz(-1.9907711) q[0];
sx q[0];
rz(-2.6429667) q[0];
rz(-pi) q[1];
rz(-2.6996428) q[2];
sx q[2];
rz(-2.0734678) q[2];
sx q[2];
rz(1.82077) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.8173302) q[1];
sx q[1];
rz(-1.8166537) q[1];
sx q[1];
rz(-0.27995963) q[1];
rz(-pi) q[2];
rz(1.1340895) q[3];
sx q[3];
rz(-2.1434382) q[3];
sx q[3];
rz(0.2084032) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.87172047) q[2];
sx q[2];
rz(-1.6951268) q[2];
sx q[2];
rz(0.63228697) q[2];
rz(-0.91442433) q[3];
sx q[3];
rz(-1.881003) q[3];
sx q[3];
rz(-2.9818592) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3226427) q[0];
sx q[0];
rz(-0.85700789) q[0];
sx q[0];
rz(-2.6981165) q[0];
rz(0.30287287) q[1];
sx q[1];
rz(-2.5587406) q[1];
sx q[1];
rz(1.3076967) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0118367) q[0];
sx q[0];
rz(-0.34809732) q[0];
sx q[0];
rz(2.6231717) q[0];
rz(-pi) q[1];
rz(0.52438122) q[2];
sx q[2];
rz(-0.57106995) q[2];
sx q[2];
rz(-1.8586161) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9759226) q[1];
sx q[1];
rz(-2.5745086) q[1];
sx q[1];
rz(2.0815064) q[1];
rz(-pi) q[2];
rz(2.0409706) q[3];
sx q[3];
rz(-1.8368145) q[3];
sx q[3];
rz(3.0373552) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.12094721) q[2];
sx q[2];
rz(-1.0409313) q[2];
sx q[2];
rz(-0.61332235) q[2];
rz(-0.81765085) q[3];
sx q[3];
rz(-2.652467) q[3];
sx q[3];
rz(0.31942719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9870975) q[0];
sx q[0];
rz(-1.7648062) q[0];
sx q[0];
rz(-2.7823271) q[0];
rz(-2.477395) q[1];
sx q[1];
rz(-1.8134873) q[1];
sx q[1];
rz(-0.42116234) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4182868) q[0];
sx q[0];
rz(-1.3483185) q[0];
sx q[0];
rz(-1.9487593) q[0];
x q[1];
rz(0.0053275544) q[2];
sx q[2];
rz(-1.498328) q[2];
sx q[2];
rz(2.7203512) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.87766552) q[1];
sx q[1];
rz(-1.8322199) q[1];
sx q[1];
rz(-1.1007453) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1499885) q[3];
sx q[3];
rz(-1.0458071) q[3];
sx q[3];
rz(2.2158509) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.26587129) q[2];
sx q[2];
rz(-3.044812) q[2];
sx q[2];
rz(1.9270012) q[2];
rz(0.38861361) q[3];
sx q[3];
rz(-0.92258421) q[3];
sx q[3];
rz(1.0987561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6471967) q[0];
sx q[0];
rz(-1.5460486) q[0];
sx q[0];
rz(-1.8478951) q[0];
rz(-2.9248059) q[1];
sx q[1];
rz(-2.4516791) q[1];
sx q[1];
rz(-2.6286415) q[1];
rz(0.51914712) q[2];
sx q[2];
rz(-0.74719528) q[2];
sx q[2];
rz(-2.0120646) q[2];
rz(-2.6445845) q[3];
sx q[3];
rz(-2.6217032) q[3];
sx q[3];
rz(1.1006127) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
