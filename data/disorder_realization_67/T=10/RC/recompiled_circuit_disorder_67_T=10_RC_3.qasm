OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.6131634) q[0];
sx q[0];
rz(-2.0818721) q[0];
sx q[0];
rz(-0.73097316) q[0];
rz(1.641474) q[1];
sx q[1];
rz(-1.0348231) q[1];
sx q[1];
rz(2.1980481) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7109414) q[0];
sx q[0];
rz(-1.5358155) q[0];
sx q[0];
rz(3.0745688) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.11124723) q[2];
sx q[2];
rz(-0.6539549) q[2];
sx q[2];
rz(0.36270579) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4409677) q[1];
sx q[1];
rz(-1.1807627) q[1];
sx q[1];
rz(2.7729211) q[1];
rz(-pi) q[2];
rz(1.7581851) q[3];
sx q[3];
rz(-1.870578) q[3];
sx q[3];
rz(-0.0084458394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.25508183) q[2];
sx q[2];
rz(-1.3811029) q[2];
sx q[2];
rz(1.8908267) q[2];
rz(1.4261774) q[3];
sx q[3];
rz(-0.91606796) q[3];
sx q[3];
rz(-2.1616518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0086867) q[0];
sx q[0];
rz(-1.0722906) q[0];
sx q[0];
rz(-0.59894484) q[0];
rz(-1.3409746) q[1];
sx q[1];
rz(-0.95021617) q[1];
sx q[1];
rz(0.96639955) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.37106284) q[0];
sx q[0];
rz(-0.75936985) q[0];
sx q[0];
rz(2.4735527) q[0];
rz(-pi) q[1];
rz(-1.6205377) q[2];
sx q[2];
rz(-0.94919862) q[2];
sx q[2];
rz(-3.0034686) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.4557138) q[1];
sx q[1];
rz(-0.13882942) q[1];
sx q[1];
rz(-2.5338737) q[1];
x q[2];
rz(-2.0004683) q[3];
sx q[3];
rz(-0.86887348) q[3];
sx q[3];
rz(-2.7817291) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.085658375) q[2];
sx q[2];
rz(-2.3264383) q[2];
sx q[2];
rz(2.8404964) q[2];
rz(-1.9484693) q[3];
sx q[3];
rz(-1.5501225) q[3];
sx q[3];
rz(1.955207) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.10467228) q[0];
sx q[0];
rz(-1.6372697) q[0];
sx q[0];
rz(-1.1874636) q[0];
rz(-1.9056412) q[1];
sx q[1];
rz(-2.1042447) q[1];
sx q[1];
rz(1.8240066) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0061958) q[0];
sx q[0];
rz(-2.4164696) q[0];
sx q[0];
rz(-0.32307415) q[0];
rz(-pi) q[1];
rz(-1.3733528) q[2];
sx q[2];
rz(-1.3609481) q[2];
sx q[2];
rz(-2.2270577) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.6685851) q[1];
sx q[1];
rz(-2.7088532) q[1];
sx q[1];
rz(0.24344484) q[1];
rz(-pi) q[2];
rz(-1.0533603) q[3];
sx q[3];
rz(-1.468717) q[3];
sx q[3];
rz(1.3641588) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1304156) q[2];
sx q[2];
rz(-1.7480363) q[2];
sx q[2];
rz(-2.9349566) q[2];
rz(-2.4335499) q[3];
sx q[3];
rz(-2.9338624) q[3];
sx q[3];
rz(-0.89282435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3669423) q[0];
sx q[0];
rz(-0.21454021) q[0];
sx q[0];
rz(2.2553717) q[0];
rz(-2.1318502) q[1];
sx q[1];
rz(-2.2354398) q[1];
sx q[1];
rz(-1.2264235) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3564295) q[0];
sx q[0];
rz(-0.7298846) q[0];
sx q[0];
rz(1.6701783) q[0];
rz(-1.0272155) q[2];
sx q[2];
rz(-0.76374861) q[2];
sx q[2];
rz(-0.56994146) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.9981421) q[1];
sx q[1];
rz(-1.3610024) q[1];
sx q[1];
rz(-1.3694805) q[1];
x q[2];
rz(-0.48638101) q[3];
sx q[3];
rz(-1.2380621) q[3];
sx q[3];
rz(-0.90852028) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.6440789) q[2];
sx q[2];
rz(-1.3138694) q[2];
sx q[2];
rz(-0.99299661) q[2];
rz(-1.8289061) q[3];
sx q[3];
rz(-1.0220746) q[3];
sx q[3];
rz(-0.48373568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.118367) q[0];
sx q[0];
rz(-2.826773) q[0];
sx q[0];
rz(1.9556048) q[0];
rz(1.7182619) q[1];
sx q[1];
rz(-1.6512197) q[1];
sx q[1];
rz(0.58247724) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5486149) q[0];
sx q[0];
rz(-0.60615221) q[0];
sx q[0];
rz(1.6968326) q[0];
x q[1];
rz(-1.3864473) q[2];
sx q[2];
rz(-0.68612387) q[2];
sx q[2];
rz(2.9550936) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9936258) q[1];
sx q[1];
rz(-2.2851351) q[1];
sx q[1];
rz(1.5961958) q[1];
rz(-2.7353103) q[3];
sx q[3];
rz(-1.5092106) q[3];
sx q[3];
rz(-0.53698925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.65486583) q[2];
sx q[2];
rz(-1.606769) q[2];
sx q[2];
rz(-0.081710903) q[2];
rz(-0.47406667) q[3];
sx q[3];
rz(-1.877955) q[3];
sx q[3];
rz(-1.4985532) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.896647) q[0];
sx q[0];
rz(-1.1279673) q[0];
sx q[0];
rz(0.0078049302) q[0];
rz(1.4004978) q[1];
sx q[1];
rz(-0.85406071) q[1];
sx q[1];
rz(2.0369464) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1305599) q[0];
sx q[0];
rz(-1.4078119) q[0];
sx q[0];
rz(-2.4736604) q[0];
rz(-pi) q[1];
rz(-2.5014624) q[2];
sx q[2];
rz(-1.374561) q[2];
sx q[2];
rz(0.093984691) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0749952) q[1];
sx q[1];
rz(-2.1146333) q[1];
sx q[1];
rz(0.72504136) q[1];
rz(-2.036318) q[3];
sx q[3];
rz(-1.421531) q[3];
sx q[3];
rz(2.535459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2509987) q[2];
sx q[2];
rz(-2.4119174) q[2];
sx q[2];
rz(0.55523038) q[2];
rz(-0.17272078) q[3];
sx q[3];
rz(-1.3135066) q[3];
sx q[3];
rz(-1.5482607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5531439) q[0];
sx q[0];
rz(-0.48184904) q[0];
sx q[0];
rz(3.0798262) q[0];
rz(0.24208367) q[1];
sx q[1];
rz(-2.7663019) q[1];
sx q[1];
rz(2.0297091) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0430849) q[0];
sx q[0];
rz(-0.83139172) q[0];
sx q[0];
rz(0.99494536) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3948453) q[2];
sx q[2];
rz(-0.07882747) q[2];
sx q[2];
rz(-1.4836756) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.4278533) q[1];
sx q[1];
rz(-2.3728328) q[1];
sx q[1];
rz(0.40366918) q[1];
rz(-pi) q[2];
rz(-1.816733) q[3];
sx q[3];
rz(-1.2517559) q[3];
sx q[3];
rz(0.46003534) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.8093439) q[2];
sx q[2];
rz(-0.76433864) q[2];
sx q[2];
rz(-2.3279482) q[2];
rz(1.7371197) q[3];
sx q[3];
rz(-0.22870326) q[3];
sx q[3];
rz(2.5261734) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(0.62548816) q[0];
sx q[0];
rz(-1.7913211) q[0];
sx q[0];
rz(1.7161436) q[0];
rz(1.5215993) q[1];
sx q[1];
rz(-0.75192538) q[1];
sx q[1];
rz(-0.63751784) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.051415074) q[0];
sx q[0];
rz(-2.1438103) q[0];
sx q[0];
rz(2.2374002) q[0];
rz(-0.39890639) q[2];
sx q[2];
rz(-2.3762694) q[2];
sx q[2];
rz(-1.4590291) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.058640826) q[1];
sx q[1];
rz(-2.4070027) q[1];
sx q[1];
rz(2.525108) q[1];
rz(-pi) q[2];
rz(2.5708837) q[3];
sx q[3];
rz(-1.686704) q[3];
sx q[3];
rz(-1.7403719) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0104388) q[2];
sx q[2];
rz(-0.70019478) q[2];
sx q[2];
rz(-1.1361702) q[2];
rz(1.4853959) q[3];
sx q[3];
rz(-0.52432004) q[3];
sx q[3];
rz(-1.9320528) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6256325) q[0];
sx q[0];
rz(-0.78173286) q[0];
sx q[0];
rz(0.67614722) q[0];
rz(2.3162084) q[1];
sx q[1];
rz(-2.717658) q[1];
sx q[1];
rz(1.0151781) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0494173) q[0];
sx q[0];
rz(-1.9855238) q[0];
sx q[0];
rz(1.7134922) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9416111) q[2];
sx q[2];
rz(-2.4790384) q[2];
sx q[2];
rz(2.8560864) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.79267348) q[1];
sx q[1];
rz(-1.9209849) q[1];
sx q[1];
rz(-2.7677571) q[1];
rz(-pi) q[2];
rz(1.2328524) q[3];
sx q[3];
rz(-2.3082808) q[3];
sx q[3];
rz(-1.0271219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4201346) q[2];
sx q[2];
rz(-0.2162424) q[2];
sx q[2];
rz(-0.96735111) q[2];
rz(1.597065) q[3];
sx q[3];
rz(-1.3128076) q[3];
sx q[3];
rz(2.7887204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-0.6185146) q[0];
sx q[0];
rz(-1.5788364) q[0];
sx q[0];
rz(0.05649795) q[0];
rz(-1.0724732) q[1];
sx q[1];
rz(-1.871855) q[1];
sx q[1];
rz(-1.7369695) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0987941) q[0];
sx q[0];
rz(-0.3846752) q[0];
sx q[0];
rz(-1.2205475) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8374412) q[2];
sx q[2];
rz(-1.6870105) q[2];
sx q[2];
rz(0.2955557) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.2056634) q[1];
sx q[1];
rz(-0.91055369) q[1];
sx q[1];
rz(-0.6026938) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0839755) q[3];
sx q[3];
rz(-2.9908097) q[3];
sx q[3];
rz(0.57612102) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.848032) q[2];
sx q[2];
rz(-2.3876987) q[2];
sx q[2];
rz(0.30612293) q[2];
rz(-0.39811578) q[3];
sx q[3];
rz(-1.6653776) q[3];
sx q[3];
rz(-1.0242296) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8067779) q[0];
sx q[0];
rz(-1.0659185) q[0];
sx q[0];
rz(0.48068) q[0];
rz(1.9819992) q[1];
sx q[1];
rz(-1.0212785) q[1];
sx q[1];
rz(-2.8550128) q[1];
rz(1.8691312) q[2];
sx q[2];
rz(-1.8229501) q[2];
sx q[2];
rz(2.0291871) q[2];
rz(-2.1383193) q[3];
sx q[3];
rz(-2.5732187) q[3];
sx q[3];
rz(-0.10789286) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
