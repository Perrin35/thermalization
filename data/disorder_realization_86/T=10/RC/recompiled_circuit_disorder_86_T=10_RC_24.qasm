OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(1.7744301) q[0];
sx q[0];
rz(-0.91355938) q[0];
sx q[0];
rz(1.4120742) q[0];
rz(-2.9867759) q[1];
sx q[1];
rz(5.6875416) q[1];
sx q[1];
rz(7.7653801) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94443653) q[0];
sx q[0];
rz(-1.7261788) q[0];
sx q[0];
rz(-2.7128501) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8843845) q[2];
sx q[2];
rz(-1.4423941) q[2];
sx q[2];
rz(-2.6103225) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.23956242) q[1];
sx q[1];
rz(-2.4671577) q[1];
sx q[1];
rz(-0.85427888) q[1];
rz(-pi) q[2];
rz(-1.7300624) q[3];
sx q[3];
rz(-2.5730238) q[3];
sx q[3];
rz(2.0025314) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.1564864) q[2];
sx q[2];
rz(-0.50922314) q[2];
sx q[2];
rz(-2.2757754) q[2];
rz(0.95430294) q[3];
sx q[3];
rz(-1.6031957) q[3];
sx q[3];
rz(1.2877864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.99825478) q[0];
sx q[0];
rz(-1.4366432) q[0];
sx q[0];
rz(-3.1153733) q[0];
rz(-1.6014618) q[1];
sx q[1];
rz(-1.5427579) q[1];
sx q[1];
rz(-2.1781133) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.547077) q[0];
sx q[0];
rz(-1.5684959) q[0];
sx q[0];
rz(-0.95675795) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.0785525) q[2];
sx q[2];
rz(-0.62765593) q[2];
sx q[2];
rz(-0.17222675) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7746437) q[1];
sx q[1];
rz(-2.0683214) q[1];
sx q[1];
rz(-2.7760387) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.96879949) q[3];
sx q[3];
rz(-2.7235944) q[3];
sx q[3];
rz(-0.19817643) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.5144689) q[2];
sx q[2];
rz(-1.1273948) q[2];
sx q[2];
rz(-3.0070686) q[2];
rz(0.7450122) q[3];
sx q[3];
rz(-0.22694215) q[3];
sx q[3];
rz(-2.1988595) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9298252) q[0];
sx q[0];
rz(-2.7524502) q[0];
sx q[0];
rz(2.3441558) q[0];
rz(-1.047661) q[1];
sx q[1];
rz(-2.9918549) q[1];
sx q[1];
rz(0.55999666) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13941923) q[0];
sx q[0];
rz(-1.1886485) q[0];
sx q[0];
rz(0.21811534) q[0];
rz(-pi) q[1];
rz(-3.0736018) q[2];
sx q[2];
rz(-0.30390856) q[2];
sx q[2];
rz(1.7169203) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.8289889) q[1];
sx q[1];
rz(-1.3723515) q[1];
sx q[1];
rz(2.7752084) q[1];
x q[2];
rz(-1.0873763) q[3];
sx q[3];
rz(-1.9994352) q[3];
sx q[3];
rz(0.40070686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.75227633) q[2];
sx q[2];
rz(-1.1976778) q[2];
sx q[2];
rz(-0.17253549) q[2];
rz(0.98207384) q[3];
sx q[3];
rz(-1.3970102) q[3];
sx q[3];
rz(-2.0836232) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1383706) q[0];
sx q[0];
rz(-2.0439742) q[0];
sx q[0];
rz(-2.8570783) q[0];
rz(0.31670397) q[1];
sx q[1];
rz(-0.43276325) q[1];
sx q[1];
rz(-1.8428615) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38628681) q[0];
sx q[0];
rz(-0.55314976) q[0];
sx q[0];
rz(-2.0095216) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.1403055) q[2];
sx q[2];
rz(-0.80438559) q[2];
sx q[2];
rz(-0.019891642) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.1620996) q[1];
sx q[1];
rz(-2.5433308) q[1];
sx q[1];
rz(1.23566) q[1];
rz(-pi) q[2];
rz(3.0292547) q[3];
sx q[3];
rz(-1.2445407) q[3];
sx q[3];
rz(-0.60409594) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5359042) q[2];
sx q[2];
rz(-2.9457592) q[2];
sx q[2];
rz(2.7569125) q[2];
rz(-0.7540594) q[3];
sx q[3];
rz(-1.056517) q[3];
sx q[3];
rz(-1.6872905) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6032747) q[0];
sx q[0];
rz(-2.1544927) q[0];
sx q[0];
rz(-1.7549365) q[0];
rz(-0.23100135) q[1];
sx q[1];
rz(-1.8004386) q[1];
sx q[1];
rz(0.2968266) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7291527) q[0];
sx q[0];
rz(-1.531633) q[0];
sx q[0];
rz(0.57106437) q[0];
rz(-pi) q[1];
rz(1.7237687) q[2];
sx q[2];
rz(-0.83380552) q[2];
sx q[2];
rz(2.134915) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.27767147) q[1];
sx q[1];
rz(-1.1754981) q[1];
sx q[1];
rz(1.0635832) q[1];
rz(-pi) q[2];
x q[2];
rz(-3.0389901) q[3];
sx q[3];
rz(-2.4253064) q[3];
sx q[3];
rz(-2.7308381) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.37830535) q[2];
sx q[2];
rz(-1.3102691) q[2];
sx q[2];
rz(-2.7491167) q[2];
rz(-1.9893507) q[3];
sx q[3];
rz(-2.4270054) q[3];
sx q[3];
rz(-0.31744441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1257989) q[0];
sx q[0];
rz(-1.5725461) q[0];
sx q[0];
rz(-0.75138599) q[0];
rz(-1.8136576) q[1];
sx q[1];
rz(-1.8782047) q[1];
sx q[1];
rz(-0.60633916) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7327001) q[0];
sx q[0];
rz(-3.0928287) q[0];
sx q[0];
rz(2.7932037) q[0];
x q[1];
rz(2.2491127) q[2];
sx q[2];
rz(-1.2391029) q[2];
sx q[2];
rz(-1.4075556) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.11583466) q[1];
sx q[1];
rz(-2.5587213) q[1];
sx q[1];
rz(1.9028266) q[1];
x q[2];
rz(1.353225) q[3];
sx q[3];
rz(-1.6440344) q[3];
sx q[3];
rz(0.53461246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.5027344) q[2];
sx q[2];
rz(-2.0998462) q[2];
sx q[2];
rz(1.139337) q[2];
rz(-1.6566488) q[3];
sx q[3];
rz(-1.1805725) q[3];
sx q[3];
rz(3.0373354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.5320324) q[0];
sx q[0];
rz(-2.39344) q[0];
sx q[0];
rz(-2.6334921) q[0];
rz(1.5628901) q[1];
sx q[1];
rz(-1.0888313) q[1];
sx q[1];
rz(2.3513444) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9433141) q[0];
sx q[0];
rz(-0.6753079) q[0];
sx q[0];
rz(-2.6576463) q[0];
rz(3.087567) q[2];
sx q[2];
rz(-0.21394357) q[2];
sx q[2];
rz(-0.33188785) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8040647) q[1];
sx q[1];
rz(-1.3009326) q[1];
sx q[1];
rz(0.40562628) q[1];
x q[2];
rz(2.4942057) q[3];
sx q[3];
rz(-1.8772519) q[3];
sx q[3];
rz(0.90741531) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.053085176) q[2];
sx q[2];
rz(-0.44184703) q[2];
sx q[2];
rz(1.4132168) q[2];
rz(-1.6648071) q[3];
sx q[3];
rz(-1.0352742) q[3];
sx q[3];
rz(-3.0800381) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4247894) q[0];
sx q[0];
rz(-3.1112818) q[0];
sx q[0];
rz(2.0943663) q[0];
rz(0.60910243) q[1];
sx q[1];
rz(-1.4139688) q[1];
sx q[1];
rz(1.75288) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7097276) q[0];
sx q[0];
rz(-2.1301529) q[0];
sx q[0];
rz(-3.1064242) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.6128204) q[2];
sx q[2];
rz(-1.1985589) q[2];
sx q[2];
rz(-2.1793274) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.82909225) q[1];
sx q[1];
rz(-1.6626076) q[1];
sx q[1];
rz(1.5593668) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1771168) q[3];
sx q[3];
rz(-1.8165605) q[3];
sx q[3];
rz(-0.59734674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.9528815) q[2];
sx q[2];
rz(-0.41023508) q[2];
sx q[2];
rz(0.88225538) q[2];
rz(1.4011718) q[3];
sx q[3];
rz(-1.1663576) q[3];
sx q[3];
rz(-1.2004948) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8700478) q[0];
sx q[0];
rz(-0.4168059) q[0];
sx q[0];
rz(-1.4260938) q[0];
rz(3.0601314) q[1];
sx q[1];
rz(-1.9790244) q[1];
sx q[1];
rz(-2.5833599) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4687913) q[0];
sx q[0];
rz(-1.1835915) q[0];
sx q[0];
rz(-2.0331435) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0874861) q[2];
sx q[2];
rz(-1.4634842) q[2];
sx q[2];
rz(-0.12003128) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.12802943) q[1];
sx q[1];
rz(-2.3131144) q[1];
sx q[1];
rz(0.086704266) q[1];
rz(-pi) q[2];
x q[2];
rz(2.6551412) q[3];
sx q[3];
rz(-1.4064186) q[3];
sx q[3];
rz(-2.1742976) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2925064) q[2];
sx q[2];
rz(-1.8722653) q[2];
sx q[2];
rz(2.9294087) q[2];
rz(2.9296181) q[3];
sx q[3];
rz(-0.68325716) q[3];
sx q[3];
rz(-1.2020483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.50487173) q[0];
sx q[0];
rz(-2.3175406) q[0];
sx q[0];
rz(1.6037534) q[0];
rz(2.3161855) q[1];
sx q[1];
rz(-2.4688265) q[1];
sx q[1];
rz(-0.5232946) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.035318035) q[0];
sx q[0];
rz(-0.9629074) q[0];
sx q[0];
rz(1.4950698) q[0];
rz(-0.043838219) q[2];
sx q[2];
rz(-2.1423116) q[2];
sx q[2];
rz(0.29585719) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.1112422) q[1];
sx q[1];
rz(-2.6551464) q[1];
sx q[1];
rz(-0.72279795) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.75533112) q[3];
sx q[3];
rz(-0.33077251) q[3];
sx q[3];
rz(-1.8473234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.3045197) q[2];
sx q[2];
rz(-0.7154811) q[2];
sx q[2];
rz(-0.26930299) q[2];
rz(0.4942016) q[3];
sx q[3];
rz(-2.2952357) q[3];
sx q[3];
rz(1.0860898) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8158648) q[0];
sx q[0];
rz(-1.6115191) q[0];
sx q[0];
rz(1.4900526) q[0];
rz(1.6745463) q[1];
sx q[1];
rz(-0.29232262) q[1];
sx q[1];
rz(1.2437337) q[1];
rz(-1.5291443) q[2];
sx q[2];
rz(-0.26298444) q[2];
sx q[2];
rz(2.0900805) q[2];
rz(1.3310824) q[3];
sx q[3];
rz(-0.78835434) q[3];
sx q[3];
rz(2.2314856) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];