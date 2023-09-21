OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.22566158) q[0];
sx q[0];
rz(-2.2731279) q[0];
sx q[0];
rz(-2.948569) q[0];
rz(1.141619) q[1];
sx q[1];
rz(-0.42998278) q[1];
sx q[1];
rz(-0.68312445) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2642759) q[0];
sx q[0];
rz(-1.6289992) q[0];
sx q[0];
rz(3.0741865) q[0];
x q[1];
rz(-1.921466) q[2];
sx q[2];
rz(-1.7979243) q[2];
sx q[2];
rz(0.35958689) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5176757) q[1];
sx q[1];
rz(-0.78849925) q[1];
sx q[1];
rz(-0.6448402) q[1];
x q[2];
rz(-2.1360374) q[3];
sx q[3];
rz(-1.894265) q[3];
sx q[3];
rz(2.2068057) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0156988) q[2];
sx q[2];
rz(-1.3796207) q[2];
sx q[2];
rz(-2.0430298) q[2];
rz(1.0788318) q[3];
sx q[3];
rz(-2.1964985) q[3];
sx q[3];
rz(2.0400955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.1612448) q[0];
sx q[0];
rz(-0.315936) q[0];
sx q[0];
rz(2.9336477) q[0];
rz(-2.5646599) q[1];
sx q[1];
rz(-2.2536342) q[1];
sx q[1];
rz(-1.6764486) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6786878) q[0];
sx q[0];
rz(-1.573274) q[0];
sx q[0];
rz(-0.72420995) q[0];
rz(-pi) q[1];
x q[1];
rz(2.8855578) q[2];
sx q[2];
rz(-1.0962152) q[2];
sx q[2];
rz(-1.408996) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.84491731) q[1];
sx q[1];
rz(-2.1293318) q[1];
sx q[1];
rz(-1.0122453) q[1];
rz(-pi) q[2];
x q[2];
rz(0.37104718) q[3];
sx q[3];
rz(-2.3608532) q[3];
sx q[3];
rz(-2.3649529) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0186105) q[2];
sx q[2];
rz(-2.1554422) q[2];
sx q[2];
rz(0.1097651) q[2];
rz(-2.5189853) q[3];
sx q[3];
rz(-0.37125769) q[3];
sx q[3];
rz(-1.4573147) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1784172) q[0];
sx q[0];
rz(-2.2816179) q[0];
sx q[0];
rz(-0.51613581) q[0];
rz(-0.57488817) q[1];
sx q[1];
rz(-2.2153885) q[1];
sx q[1];
rz(-0.80054545) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1313989) q[0];
sx q[0];
rz(-0.43361615) q[0];
sx q[0];
rz(-1.2244768) q[0];
rz(-pi) q[1];
rz(-0.86513743) q[2];
sx q[2];
rz(-1.0312928) q[2];
sx q[2];
rz(-1.3441966) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.264818) q[1];
sx q[1];
rz(-0.41300981) q[1];
sx q[1];
rz(3.1289711) q[1];
rz(-pi) q[2];
rz(-0.3018467) q[3];
sx q[3];
rz(-0.89248025) q[3];
sx q[3];
rz(0.18920004) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4642554) q[2];
sx q[2];
rz(-2.8223473) q[2];
sx q[2];
rz(1.3595954) q[2];
rz(-2.1740186) q[3];
sx q[3];
rz(-1.273497) q[3];
sx q[3];
rz(1.7165855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1490705) q[0];
sx q[0];
rz(-1.8795805) q[0];
sx q[0];
rz(-0.46491369) q[0];
rz(0.34856302) q[1];
sx q[1];
rz(-0.26270738) q[1];
sx q[1];
rz(1.0850614) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5269055) q[0];
sx q[0];
rz(-0.58262107) q[0];
sx q[0];
rz(2.7132062) q[0];
x q[1];
rz(-2.1539139) q[2];
sx q[2];
rz(-2.7104125) q[2];
sx q[2];
rz(1.093986) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.3738721) q[1];
sx q[1];
rz(-1.4704736) q[1];
sx q[1];
rz(-0.14212455) q[1];
x q[2];
rz(-0.61120175) q[3];
sx q[3];
rz(-1.2308321) q[3];
sx q[3];
rz(1.1838278) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.00502914) q[2];
sx q[2];
rz(-2.0932784) q[2];
sx q[2];
rz(-2.4604649) q[2];
rz(-2.629771) q[3];
sx q[3];
rz(-2.8183283) q[3];
sx q[3];
rz(-2.8779023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
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
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8624449) q[0];
sx q[0];
rz(-0.537916) q[0];
sx q[0];
rz(-1.408668) q[0];
rz(0.43235835) q[1];
sx q[1];
rz(-2.2996348) q[1];
sx q[1];
rz(2.1599105) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0911134) q[0];
sx q[0];
rz(-1.0266773) q[0];
sx q[0];
rz(1.0247466) q[0];
rz(-pi) q[1];
rz(1.328674) q[2];
sx q[2];
rz(-1.4171346) q[2];
sx q[2];
rz(-0.15649934) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.423646) q[1];
sx q[1];
rz(-1.4232993) q[1];
sx q[1];
rz(0.46848483) q[1];
rz(-pi) q[2];
rz(-1.4622299) q[3];
sx q[3];
rz(-2.541399) q[3];
sx q[3];
rz(-1.1443646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.1061873) q[2];
sx q[2];
rz(-1.4865439) q[2];
sx q[2];
rz(0.48941082) q[2];
rz(2.126746) q[3];
sx q[3];
rz(-2.615052) q[3];
sx q[3];
rz(1.3283407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6569825) q[0];
sx q[0];
rz(-2.9841612) q[0];
sx q[0];
rz(-2.7691675) q[0];
rz(1.3308446) q[1];
sx q[1];
rz(-1.0354038) q[1];
sx q[1];
rz(2.9763124) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3344791) q[0];
sx q[0];
rz(-1.4757336) q[0];
sx q[0];
rz(-1.6001742) q[0];
rz(-pi) q[1];
rz(-0.084947305) q[2];
sx q[2];
rz(-2.4521378) q[2];
sx q[2];
rz(1.8817608) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2416934) q[1];
sx q[1];
rz(-0.69679931) q[1];
sx q[1];
rz(-2.4844869) q[1];
x q[2];
rz(2.0323456) q[3];
sx q[3];
rz(-0.49406067) q[3];
sx q[3];
rz(0.1915313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.91339397) q[2];
sx q[2];
rz(-1.0281576) q[2];
sx q[2];
rz(1.6983263) q[2];
rz(-2.7741487) q[3];
sx q[3];
rz(-1.2747217) q[3];
sx q[3];
rz(2.929556) q[3];
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
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7845602) q[0];
sx q[0];
rz(-2.050188) q[0];
sx q[0];
rz(2.7923287) q[0];
rz(-0.7473942) q[1];
sx q[1];
rz(-2.8458197) q[1];
sx q[1];
rz(-2.4051037) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.38024494) q[0];
sx q[0];
rz(-1.5537098) q[0];
sx q[0];
rz(1.5218309) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9095124) q[2];
sx q[2];
rz(-1.3352397) q[2];
sx q[2];
rz(2.2155857) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.4941102) q[1];
sx q[1];
rz(-2.3586015) q[1];
sx q[1];
rz(-0.25581911) q[1];
rz(0.50076671) q[3];
sx q[3];
rz(-1.9541249) q[3];
sx q[3];
rz(-0.70481833) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0050469) q[2];
sx q[2];
rz(-1.9133291) q[2];
sx q[2];
rz(0.53156701) q[2];
rz(-0.68938869) q[3];
sx q[3];
rz(-1.6770984) q[3];
sx q[3];
rz(-1.846107) q[3];
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
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.078995973) q[0];
sx q[0];
rz(-1.2051219) q[0];
sx q[0];
rz(1.8435562) q[0];
rz(-2.334306) q[1];
sx q[1];
rz(-1.9629982) q[1];
sx q[1];
rz(0.92179006) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4661134) q[0];
sx q[0];
rz(-2.4371494) q[0];
sx q[0];
rz(0.2091614) q[0];
rz(-pi) q[1];
rz(2.8554101) q[2];
sx q[2];
rz(-0.85368644) q[2];
sx q[2];
rz(-1.66695) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.39107716) q[1];
sx q[1];
rz(-0.32954307) q[1];
sx q[1];
rz(-2.4604172) q[1];
rz(0.17955762) q[3];
sx q[3];
rz(-0.87053821) q[3];
sx q[3];
rz(1.8732656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.2395997) q[2];
sx q[2];
rz(-2.5805876) q[2];
sx q[2];
rz(-1.1716589) q[2];
rz(1.7840067) q[3];
sx q[3];
rz(-1.4566908) q[3];
sx q[3];
rz(1.4484891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8885324) q[0];
sx q[0];
rz(-1.9655515) q[0];
sx q[0];
rz(1.4755479) q[0];
rz(-1.6015923) q[1];
sx q[1];
rz(-1.4614636) q[1];
sx q[1];
rz(-2.4618861) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8863227) q[0];
sx q[0];
rz(-2.5810044) q[0];
sx q[0];
rz(-0.38829304) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.26608) q[2];
sx q[2];
rz(-1.4317703) q[2];
sx q[2];
rz(-1.306844) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.3956086) q[1];
sx q[1];
rz(-1.9095699) q[1];
sx q[1];
rz(2.4376274) q[1];
rz(-pi) q[2];
rz(-2.4623975) q[3];
sx q[3];
rz(-1.067357) q[3];
sx q[3];
rz(-0.60590832) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1264964) q[2];
sx q[2];
rz(-1.482778) q[2];
sx q[2];
rz(2.2311907) q[2];
rz(-0.67534584) q[3];
sx q[3];
rz(-0.93674913) q[3];
sx q[3];
rz(-2.2909686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.05242059) q[0];
sx q[0];
rz(-1.9071254) q[0];
sx q[0];
rz(-2.4269379) q[0];
rz(2.4275298) q[1];
sx q[1];
rz(-0.95497447) q[1];
sx q[1];
rz(1.1766599) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6274174) q[0];
sx q[0];
rz(-1.5222933) q[0];
sx q[0];
rz(1.4072627) q[0];
rz(-pi) q[1];
x q[1];
rz(0.87877019) q[2];
sx q[2];
rz(-0.76239097) q[2];
sx q[2];
rz(-2.9825485) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9709819) q[1];
sx q[1];
rz(-2.5787528) q[1];
sx q[1];
rz(1.4303722) q[1];
rz(1.482974) q[3];
sx q[3];
rz(-2.0801968) q[3];
sx q[3];
rz(2.1109964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8979793) q[2];
sx q[2];
rz(-1.672013) q[2];
sx q[2];
rz(-1.127355) q[2];
rz(0.35774287) q[3];
sx q[3];
rz(-0.8711516) q[3];
sx q[3];
rz(-1.0424785) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7174299) q[0];
sx q[0];
rz(-1.9308199) q[0];
sx q[0];
rz(0.45817026) q[0];
rz(-0.39623109) q[1];
sx q[1];
rz(-3.1165262) q[1];
sx q[1];
rz(-2.810626) q[1];
rz(-2.836543) q[2];
sx q[2];
rz(-2.3813644) q[2];
sx q[2];
rz(2.7347953) q[2];
rz(1.5626004) q[3];
sx q[3];
rz(-1.2200439) q[3];
sx q[3];
rz(-0.72190819) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];