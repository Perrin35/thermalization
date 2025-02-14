OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.5105628) q[0];
sx q[0];
rz(-0.50463843) q[0];
sx q[0];
rz(2.7161427) q[0];
rz(-2.5692441) q[1];
sx q[1];
rz(-1.0971789) q[1];
sx q[1];
rz(-1.2001349) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6365189) q[0];
sx q[0];
rz(-0.25101343) q[0];
sx q[0];
rz(0.039029718) q[0];
rz(2.2283353) q[2];
sx q[2];
rz(-2.7139671) q[2];
sx q[2];
rz(-2.0562125) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.8737826) q[1];
sx q[1];
rz(-1.7846196) q[1];
sx q[1];
rz(-0.38114433) q[1];
rz(-0.17961924) q[3];
sx q[3];
rz(-0.67905513) q[3];
sx q[3];
rz(-1.2509026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.1537062) q[2];
sx q[2];
rz(-1.989863) q[2];
sx q[2];
rz(-0.64725867) q[2];
rz(2.9636532) q[3];
sx q[3];
rz(-1.8837594) q[3];
sx q[3];
rz(1.2037163) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7385638) q[0];
sx q[0];
rz(-3.0572427) q[0];
sx q[0];
rz(-2.6179598) q[0];
rz(-1.9117484) q[1];
sx q[1];
rz(-2.3975027) q[1];
sx q[1];
rz(-2.8935166) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1551127) q[0];
sx q[0];
rz(-2.393828) q[0];
sx q[0];
rz(2.98481) q[0];
x q[1];
rz(2.352574) q[2];
sx q[2];
rz(-0.58525267) q[2];
sx q[2];
rz(-1.8155797) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.7404799) q[1];
sx q[1];
rz(-1.4251627) q[1];
sx q[1];
rz(3.1023953) q[1];
x q[2];
rz(2.5026191) q[3];
sx q[3];
rz(-2.3666214) q[3];
sx q[3];
rz(2.4955622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.5292042) q[2];
sx q[2];
rz(-2.2233621) q[2];
sx q[2];
rz(-0.53981346) q[2];
rz(2.9811033) q[3];
sx q[3];
rz(-1.548998) q[3];
sx q[3];
rz(2.3314355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6432583) q[0];
sx q[0];
rz(-0.67688268) q[0];
sx q[0];
rz(0.13370378) q[0];
rz(1.5191822) q[1];
sx q[1];
rz(-1.2956053) q[1];
sx q[1];
rz(-0.79230961) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3334219) q[0];
sx q[0];
rz(-1.9229888) q[0];
sx q[0];
rz(-2.9461224) q[0];
x q[1];
rz(0.14459212) q[2];
sx q[2];
rz(-1.4404313) q[2];
sx q[2];
rz(-3.0693288) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.50666675) q[1];
sx q[1];
rz(-0.95639765) q[1];
sx q[1];
rz(2.3022815) q[1];
x q[2];
rz(1.5706704) q[3];
sx q[3];
rz(-1.8709261) q[3];
sx q[3];
rz(-0.89849328) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.2866659) q[2];
sx q[2];
rz(-2.4714405) q[2];
sx q[2];
rz(-2.3112042) q[2];
rz(-1.7928127) q[3];
sx q[3];
rz(-1.034779) q[3];
sx q[3];
rz(-0.59876284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8778359) q[0];
sx q[0];
rz(-2.1649375) q[0];
sx q[0];
rz(-0.42309716) q[0];
rz(1.1571723) q[1];
sx q[1];
rz(-1.4856228) q[1];
sx q[1];
rz(0.62209904) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7579278) q[0];
sx q[0];
rz(-0.94449857) q[0];
sx q[0];
rz(-0.48547283) q[0];
rz(1.1272264) q[2];
sx q[2];
rz(-2.6831919) q[2];
sx q[2];
rz(-1.0215525) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.9378949) q[1];
sx q[1];
rz(-2.7337498) q[1];
sx q[1];
rz(-1.6126584) q[1];
rz(1.6016209) q[3];
sx q[3];
rz(-2.4158231) q[3];
sx q[3];
rz(2.9963065) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.28079924) q[2];
sx q[2];
rz(-1.7223027) q[2];
sx q[2];
rz(0.61817509) q[2];
rz(-1.482796) q[3];
sx q[3];
rz(-1.7890472) q[3];
sx q[3];
rz(-1.6331204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16267714) q[0];
sx q[0];
rz(-1.8122883) q[0];
sx q[0];
rz(2.3090114) q[0];
rz(2.816448) q[1];
sx q[1];
rz(-2.145576) q[1];
sx q[1];
rz(-0.064373374) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8125809) q[0];
sx q[0];
rz(-1.4509177) q[0];
sx q[0];
rz(-0.32696149) q[0];
rz(-pi) q[1];
rz(-2.5413248) q[2];
sx q[2];
rz(-1.9429328) q[2];
sx q[2];
rz(-2.0780217) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.87582237) q[1];
sx q[1];
rz(-1.0604621) q[1];
sx q[1];
rz(0.93632621) q[1];
x q[2];
rz(-0.35803087) q[3];
sx q[3];
rz(-2.1674986) q[3];
sx q[3];
rz(2.1880414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.23318204) q[2];
sx q[2];
rz(-0.57047129) q[2];
sx q[2];
rz(-1.2916267) q[2];
rz(-2.6455247) q[3];
sx q[3];
rz(-0.78947624) q[3];
sx q[3];
rz(1.1424278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18601501) q[0];
sx q[0];
rz(-2.8514974) q[0];
sx q[0];
rz(0.41906038) q[0];
rz(0.31173197) q[1];
sx q[1];
rz(-1.5429976) q[1];
sx q[1];
rz(2.9980803) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1868033) q[0];
sx q[0];
rz(-2.020103) q[0];
sx q[0];
rz(1.1727929) q[0];
rz(2.4359598) q[2];
sx q[2];
rz(-2.7876283) q[2];
sx q[2];
rz(1.3753381) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.011605) q[1];
sx q[1];
rz(-1.2530348) q[1];
sx q[1];
rz(-2.3684273) q[1];
x q[2];
rz(0.26007248) q[3];
sx q[3];
rz(-1.4315769) q[3];
sx q[3];
rz(-1.8782488) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.781337) q[2];
sx q[2];
rz(-2.2402253) q[2];
sx q[2];
rz(-2.0085013) q[2];
rz(1.1059149) q[3];
sx q[3];
rz(-1.4191041) q[3];
sx q[3];
rz(-2.6667986) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.356242) q[0];
sx q[0];
rz(-0.79371912) q[0];
sx q[0];
rz(-1.6957977) q[0];
rz(-2.4244507) q[1];
sx q[1];
rz(-1.278806) q[1];
sx q[1];
rz(-2.380127) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4091771) q[0];
sx q[0];
rz(-1.4237561) q[0];
sx q[0];
rz(2.5086286) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1091719) q[2];
sx q[2];
rz(-0.79029492) q[2];
sx q[2];
rz(0.99144713) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.95132212) q[1];
sx q[1];
rz(-0.99748625) q[1];
sx q[1];
rz(-3.0828397) q[1];
rz(-pi) q[2];
rz(2.981212) q[3];
sx q[3];
rz(-1.7564991) q[3];
sx q[3];
rz(-1.4726313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.7722499) q[2];
sx q[2];
rz(-2.6476634) q[2];
sx q[2];
rz(-2.210145) q[2];
rz(2.5233968) q[3];
sx q[3];
rz(-1.5074916) q[3];
sx q[3];
rz(-2.1759694) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44929993) q[0];
sx q[0];
rz(-1.7789142) q[0];
sx q[0];
rz(-1.3344673) q[0];
rz(-0.5131228) q[1];
sx q[1];
rz(-2.158439) q[1];
sx q[1];
rz(1.9122874) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2199101) q[0];
sx q[0];
rz(-1.1510885) q[0];
sx q[0];
rz(2.3786663) q[0];
rz(-pi) q[1];
rz(1.233835) q[2];
sx q[2];
rz(-2.8671727) q[2];
sx q[2];
rz(1.929259) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.25101623) q[1];
sx q[1];
rz(-2.9843247) q[1];
sx q[1];
rz(-0.1319488) q[1];
x q[2];
rz(-0.02259262) q[3];
sx q[3];
rz(-0.66231662) q[3];
sx q[3];
rz(2.3316771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.3076155) q[2];
sx q[2];
rz(-2.6145085) q[2];
sx q[2];
rz(2.9311467) q[2];
rz(0.065841913) q[3];
sx q[3];
rz(-0.47271553) q[3];
sx q[3];
rz(2.3426447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54302067) q[0];
sx q[0];
rz(-0.6544756) q[0];
sx q[0];
rz(-1.8684813) q[0];
rz(-1.8674564) q[1];
sx q[1];
rz(-2.4576063) q[1];
sx q[1];
rz(-0.88409105) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1696816) q[0];
sx q[0];
rz(-2.187664) q[0];
sx q[0];
rz(3.141298) q[0];
x q[1];
rz(-0.42914088) q[2];
sx q[2];
rz(-2.57518) q[2];
sx q[2];
rz(-0.31508101) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6940004) q[1];
sx q[1];
rz(-1.847637) q[1];
sx q[1];
rz(2.0029699) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.8330634) q[3];
sx q[3];
rz(-1.8179699) q[3];
sx q[3];
rz(0.17251523) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.990443) q[2];
sx q[2];
rz(-0.98439011) q[2];
sx q[2];
rz(-1.7669862) q[2];
rz(2.9511342) q[3];
sx q[3];
rz(-2.3951267) q[3];
sx q[3];
rz(-1.9577352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9780438) q[0];
sx q[0];
rz(-1.4530285) q[0];
sx q[0];
rz(1.2646041) q[0];
rz(3.0679852) q[1];
sx q[1];
rz(-1.0818447) q[1];
sx q[1];
rz(0.61990613) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4891751) q[0];
sx q[0];
rz(-0.76120725) q[0];
sx q[0];
rz(-1.3309782) q[0];
x q[1];
rz(-0.61087278) q[2];
sx q[2];
rz(-1.1430642) q[2];
sx q[2];
rz(-1.3547225) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.82295361) q[1];
sx q[1];
rz(-1.3711437) q[1];
sx q[1];
rz(2.9266788) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6750071) q[3];
sx q[3];
rz(-1.3562725) q[3];
sx q[3];
rz(-1.1669351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.8668883) q[2];
sx q[2];
rz(-2.4394749) q[2];
sx q[2];
rz(2.9956024) q[2];
rz(-0.22249666) q[3];
sx q[3];
rz(-0.22495088) q[3];
sx q[3];
rz(-0.72993025) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.20919007) q[0];
sx q[0];
rz(-1.9575735) q[0];
sx q[0];
rz(0.14458543) q[0];
rz(1.9119541) q[1];
sx q[1];
rz(-1.3508136) q[1];
sx q[1];
rz(-1.2523686) q[1];
rz(0.59496224) q[2];
sx q[2];
rz(-1.7528201) q[2];
sx q[2];
rz(2.5217944) q[2];
rz(-1.7810697) q[3];
sx q[3];
rz(-1.6023286) q[3];
sx q[3];
rz(-0.31983432) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
