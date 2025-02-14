OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.98874918) q[0];
sx q[0];
rz(-0.56892836) q[0];
sx q[0];
rz(-0.95175728) q[0];
rz(3.9495502) q[1];
sx q[1];
rz(6.1679975) q[1];
sx q[1];
rz(8.4313784) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4848413) q[0];
sx q[0];
rz(-0.97545058) q[0];
sx q[0];
rz(0.30466051) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1615218) q[2];
sx q[2];
rz(-0.23071846) q[2];
sx q[2];
rz(-1.7144817) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.3776855) q[1];
sx q[1];
rz(-1.0432668) q[1];
sx q[1];
rz(0.94859132) q[1];
rz(-2.4798054) q[3];
sx q[3];
rz(-1.6787623) q[3];
sx q[3];
rz(-1.770883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.6318165) q[2];
sx q[2];
rz(-0.6002554) q[2];
sx q[2];
rz(2.996345) q[2];
rz(-0.97233573) q[3];
sx q[3];
rz(-1.5147246) q[3];
sx q[3];
rz(1.2783031) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4165322) q[0];
sx q[0];
rz(-2.5711377) q[0];
sx q[0];
rz(2.3544627) q[0];
rz(-2.9540673) q[1];
sx q[1];
rz(-1.3860605) q[1];
sx q[1];
rz(-0.010201605) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1928471) q[0];
sx q[0];
rz(-1.5625573) q[0];
sx q[0];
rz(-1.6798937) q[0];
rz(-2.0554916) q[2];
sx q[2];
rz(-1.5124413) q[2];
sx q[2];
rz(-1.397152) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.82380159) q[1];
sx q[1];
rz(-1.4483869) q[1];
sx q[1];
rz(2.9106004) q[1];
rz(-2.3115029) q[3];
sx q[3];
rz(-1.4674868) q[3];
sx q[3];
rz(-1.938397) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-3.0804704) q[2];
sx q[2];
rz(-0.11961131) q[2];
sx q[2];
rz(-1.1629026) q[2];
rz(-1.2969147) q[3];
sx q[3];
rz(-1.4328522) q[3];
sx q[3];
rz(3.1356623) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6892683) q[0];
sx q[0];
rz(-0.56832123) q[0];
sx q[0];
rz(2.7962621) q[0];
rz(0.055222424) q[1];
sx q[1];
rz(-0.60631141) q[1];
sx q[1];
rz(2.9323554) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5674735) q[0];
sx q[0];
rz(-0.86025877) q[0];
sx q[0];
rz(-3.0379651) q[0];
rz(-pi) q[1];
rz(-0.018304869) q[2];
sx q[2];
rz(-1.4985111) q[2];
sx q[2];
rz(-0.56906453) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.71007939) q[1];
sx q[1];
rz(-2.4071879) q[1];
sx q[1];
rz(-2.9148352) q[1];
rz(0.24233992) q[3];
sx q[3];
rz(-2.3018357) q[3];
sx q[3];
rz(-3.1104607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0914586) q[2];
sx q[2];
rz(-1.4860934) q[2];
sx q[2];
rz(-2.6996108) q[2];
rz(0.4778536) q[3];
sx q[3];
rz(-1.7171532) q[3];
sx q[3];
rz(1.2598239) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1314989) q[0];
sx q[0];
rz(-0.57259125) q[0];
sx q[0];
rz(-1.9189438) q[0];
rz(-1.6652416) q[1];
sx q[1];
rz(-1.0226701) q[1];
sx q[1];
rz(-0.30278912) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0593392) q[0];
sx q[0];
rz(-2.6586464) q[0];
sx q[0];
rz(-0.63697908) q[0];
rz(-pi) q[1];
x q[1];
rz(0.035149375) q[2];
sx q[2];
rz(-1.6498696) q[2];
sx q[2];
rz(3.1052232) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0278182) q[1];
sx q[1];
rz(-0.31732761) q[1];
sx q[1];
rz(1.7185663) q[1];
rz(-pi) q[2];
x q[2];
rz(0.62742558) q[3];
sx q[3];
rz(-1.0194821) q[3];
sx q[3];
rz(2.3105636) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0629695) q[2];
sx q[2];
rz(-1.4457694) q[2];
sx q[2];
rz(1.7930188) q[2];
rz(0.013269987) q[3];
sx q[3];
rz(-2.477406) q[3];
sx q[3];
rz(1.2191023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9297946) q[0];
sx q[0];
rz(-0.12166611) q[0];
sx q[0];
rz(-1.1965055) q[0];
rz(-0.78367805) q[1];
sx q[1];
rz(-2.13089) q[1];
sx q[1];
rz(0.7799305) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0257192) q[0];
sx q[0];
rz(-0.78479973) q[0];
sx q[0];
rz(2.7682106) q[0];
x q[1];
rz(1.5151603) q[2];
sx q[2];
rz(-2.0075744) q[2];
sx q[2];
rz(1.7064777) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3253567) q[1];
sx q[1];
rz(-0.79918095) q[1];
sx q[1];
rz(-2.5926209) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2524917) q[3];
sx q[3];
rz(-1.0530143) q[3];
sx q[3];
rz(-0.095967231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.1448867) q[2];
sx q[2];
rz(-0.11652623) q[2];
sx q[2];
rz(-3.0987926) q[2];
rz(-1.5329817) q[3];
sx q[3];
rz(-1.3442842) q[3];
sx q[3];
rz(-2.2406421) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.75339371) q[0];
sx q[0];
rz(-2.2008984) q[0];
sx q[0];
rz(-0.70145506) q[0];
rz(-2.258621) q[1];
sx q[1];
rz(-1.7794304) q[1];
sx q[1];
rz(-2.0515474) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0018679) q[0];
sx q[0];
rz(-1.1764948) q[0];
sx q[0];
rz(0.35625881) q[0];
rz(2.0284611) q[2];
sx q[2];
rz(-0.90137611) q[2];
sx q[2];
rz(0.078851141) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.54451128) q[1];
sx q[1];
rz(-1.7131355) q[1];
sx q[1];
rz(-0.049171731) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.9771776) q[3];
sx q[3];
rz(-1.3344171) q[3];
sx q[3];
rz(-2.6316093) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.28114721) q[2];
sx q[2];
rz(-2.2522085) q[2];
sx q[2];
rz(1.7640198) q[2];
rz(0.095976202) q[3];
sx q[3];
rz(-2.4317661) q[3];
sx q[3];
rz(-0.53203741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1808566) q[0];
sx q[0];
rz(-2.2611698) q[0];
sx q[0];
rz(1.9218504) q[0];
rz(2.5324054) q[1];
sx q[1];
rz(-1.6222619) q[1];
sx q[1];
rz(-0.46527299) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2060836) q[0];
sx q[0];
rz(-0.92496189) q[0];
sx q[0];
rz(-0.95316621) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.96107595) q[2];
sx q[2];
rz(-1.3111918) q[2];
sx q[2];
rz(1.0566835) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.0445932) q[1];
sx q[1];
rz(-1.1247509) q[1];
sx q[1];
rz(3.0924132) q[1];
rz(1.3749397) q[3];
sx q[3];
rz(-0.80482414) q[3];
sx q[3];
rz(2.5087961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.51218) q[2];
sx q[2];
rz(-0.30567726) q[2];
sx q[2];
rz(-0.42210397) q[2];
rz(-3.1230538) q[3];
sx q[3];
rz(-1.5944578) q[3];
sx q[3];
rz(-0.6137994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2043532) q[0];
sx q[0];
rz(-2.6528907) q[0];
sx q[0];
rz(-0.82116425) q[0];
rz(0.69860727) q[1];
sx q[1];
rz(-0.76118529) q[1];
sx q[1];
rz(3.1331114) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9152731) q[0];
sx q[0];
rz(-0.92585671) q[0];
sx q[0];
rz(1.1862832) q[0];
rz(-pi) q[1];
x q[1];
rz(0.94363801) q[2];
sx q[2];
rz(-2.5211589) q[2];
sx q[2];
rz(2.6629184) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1084794) q[1];
sx q[1];
rz(-0.81711191) q[1];
sx q[1];
rz(-1.2014649) q[1];
rz(-pi) q[2];
rz(-0.04106122) q[3];
sx q[3];
rz(-2.2143737) q[3];
sx q[3];
rz(-0.82541556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8734133) q[2];
sx q[2];
rz(-0.11205967) q[2];
sx q[2];
rz(2.6123135) q[2];
rz(-1.4389634) q[3];
sx q[3];
rz(-1.8301423) q[3];
sx q[3];
rz(-2.8992991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5684123) q[0];
sx q[0];
rz(-2.4127164) q[0];
sx q[0];
rz(2.6045784) q[0];
rz(2.2611387) q[1];
sx q[1];
rz(-1.8388137) q[1];
sx q[1];
rz(0.74768487) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0949919) q[0];
sx q[0];
rz(-1.9032626) q[0];
sx q[0];
rz(0.73128043) q[0];
x q[1];
rz(-0.40966655) q[2];
sx q[2];
rz(-2.1958399) q[2];
sx q[2];
rz(-1.9517488) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.6612271) q[1];
sx q[1];
rz(-1.2496619) q[1];
sx q[1];
rz(0.25471806) q[1];
x q[2];
rz(-2.791312) q[3];
sx q[3];
rz(-0.64207388) q[3];
sx q[3];
rz(-2.9779899) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.65413862) q[2];
sx q[2];
rz(-0.70709252) q[2];
sx q[2];
rz(-2.1716165) q[2];
rz(-1.4966494) q[3];
sx q[3];
rz(-2.1275438) q[3];
sx q[3];
rz(1.0223201) q[3];
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
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8146166) q[0];
sx q[0];
rz(-1.6642445) q[0];
sx q[0];
rz(0.60272637) q[0];
rz(-0.0043491443) q[1];
sx q[1];
rz(-1.8252204) q[1];
sx q[1];
rz(2.0306921) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2839226) q[0];
sx q[0];
rz(-2.1956148) q[0];
sx q[0];
rz(-0.83832534) q[0];
rz(-pi) q[1];
x q[1];
rz(0.96221417) q[2];
sx q[2];
rz(-0.79691821) q[2];
sx q[2];
rz(-1.4685517) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.37392426) q[1];
sx q[1];
rz(-1.4490738) q[1];
sx q[1];
rz(-2.5472104) q[1];
x q[2];
rz(-2.6708665) q[3];
sx q[3];
rz(-1.2743055) q[3];
sx q[3];
rz(-0.52213269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.6415619) q[2];
sx q[2];
rz(-1.2189453) q[2];
sx q[2];
rz(-0.59263372) q[2];
rz(-2.5617013) q[3];
sx q[3];
rz(-0.65120828) q[3];
sx q[3];
rz(1.7013223) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92689571) q[0];
sx q[0];
rz(-1.0746645) q[0];
sx q[0];
rz(-0.97434531) q[0];
rz(0.018085619) q[1];
sx q[1];
rz(-2.7985202) q[1];
sx q[1];
rz(0.090029686) q[1];
rz(0.30822538) q[2];
sx q[2];
rz(-1.3267953) q[2];
sx q[2];
rz(0.78143668) q[2];
rz(2.3961801) q[3];
sx q[3];
rz(-1.1996307) q[3];
sx q[3];
rz(-1.2974056) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
