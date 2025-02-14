OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.087091669) q[0];
sx q[0];
rz(-2.2890685) q[0];
sx q[0];
rz(-0.26242119) q[0];
rz(-3.4497058) q[1];
sx q[1];
rz(0.83642712) q[1];
sx q[1];
rz(15.717419) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1605837) q[0];
sx q[0];
rz(-1.137059) q[0];
sx q[0];
rz(2.1061312) q[0];
x q[1];
rz(0.043536206) q[2];
sx q[2];
rz(-1.5015209) q[2];
sx q[2];
rz(-2.2486698) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5227185) q[1];
sx q[1];
rz(-0.98509231) q[1];
sx q[1];
rz(0.5962836) q[1];
rz(-pi) q[2];
rz(-2.6515048) q[3];
sx q[3];
rz(-2.2667829) q[3];
sx q[3];
rz(1.127632) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.0414163) q[2];
sx q[2];
rz(-2.2569423) q[2];
sx q[2];
rz(-2.6803988) q[2];
rz(0.50208107) q[3];
sx q[3];
rz(-2.3177948) q[3];
sx q[3];
rz(1.4095149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13040386) q[0];
sx q[0];
rz(-1.4606425) q[0];
sx q[0];
rz(-0.23656626) q[0];
rz(-0.81831167) q[1];
sx q[1];
rz(-1.9383483) q[1];
sx q[1];
rz(0.94258211) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6494076) q[0];
sx q[0];
rz(-1.5680917) q[0];
sx q[0];
rz(1.6016141) q[0];
rz(-pi) q[1];
rz(-1.4497527) q[2];
sx q[2];
rz(-1.6377875) q[2];
sx q[2];
rz(-3.071272) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.2494804) q[1];
sx q[1];
rz(-1.5693519) q[1];
sx q[1];
rz(-0.73489983) q[1];
rz(0.97901042) q[3];
sx q[3];
rz(-0.56603449) q[3];
sx q[3];
rz(2.4477959) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4054823) q[2];
sx q[2];
rz(-0.63850275) q[2];
sx q[2];
rz(-0.68518266) q[2];
rz(2.3403366) q[3];
sx q[3];
rz(-1.6359436) q[3];
sx q[3];
rz(-1.8906458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.7291229) q[0];
sx q[0];
rz(-2.6710489) q[0];
sx q[0];
rz(-0.98010081) q[0];
rz(-3.0209172) q[1];
sx q[1];
rz(-1.1680892) q[1];
sx q[1];
rz(1.4792222) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2421766) q[0];
sx q[0];
rz(-2.3606395) q[0];
sx q[0];
rz(2.910898) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1086039) q[2];
sx q[2];
rz(-2.0304907) q[2];
sx q[2];
rz(-0.60079702) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.190158) q[1];
sx q[1];
rz(-3.1312761) q[1];
sx q[1];
rz(-0.52811388) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5517518) q[3];
sx q[3];
rz(-1.2882659) q[3];
sx q[3];
rz(2.2577406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.4270619) q[2];
sx q[2];
rz(-1.6468628) q[2];
sx q[2];
rz(-2.96116) q[2];
rz(-0.42029542) q[3];
sx q[3];
rz(-0.97729483) q[3];
sx q[3];
rz(0.54496566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9412398) q[0];
sx q[0];
rz(-1.7914597) q[0];
sx q[0];
rz(3.0923162) q[0];
rz(2.409528) q[1];
sx q[1];
rz(-2.2923636) q[1];
sx q[1];
rz(1.7656309) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9889209) q[0];
sx q[0];
rz(-1.7537358) q[0];
sx q[0];
rz(-1.4519917) q[0];
rz(-pi) q[1];
rz(-0.21993262) q[2];
sx q[2];
rz(-0.77785197) q[2];
sx q[2];
rz(2.5663515) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.298784) q[1];
sx q[1];
rz(-2.7453303) q[1];
sx q[1];
rz(2.9664459) q[1];
rz(-pi) q[2];
rz(2.2729418) q[3];
sx q[3];
rz(-2.0929061) q[3];
sx q[3];
rz(-2.3412395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.10776082) q[2];
sx q[2];
rz(-1.7511448) q[2];
sx q[2];
rz(2.4793009) q[2];
rz(1.3307339) q[3];
sx q[3];
rz(-2.3190506) q[3];
sx q[3];
rz(-1.2983769) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0419643) q[0];
sx q[0];
rz(-2.8350416) q[0];
sx q[0];
rz(2.1339259) q[0];
rz(0.10534605) q[1];
sx q[1];
rz(-2.4844929) q[1];
sx q[1];
rz(-2.9109921) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4661575) q[0];
sx q[0];
rz(-2.648232) q[0];
sx q[0];
rz(-0.81146474) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.23955524) q[2];
sx q[2];
rz(-1.2192246) q[2];
sx q[2];
rz(-2.7186175) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.8523753) q[1];
sx q[1];
rz(-1.4061478) q[1];
sx q[1];
rz(1.6573805) q[1];
x q[2];
rz(-1.4197639) q[3];
sx q[3];
rz(-2.2103035) q[3];
sx q[3];
rz(3.1139656) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.0577804) q[2];
sx q[2];
rz(-0.96692204) q[2];
sx q[2];
rz(-0.18873611) q[2];
rz(-1.905929) q[3];
sx q[3];
rz(-0.50714791) q[3];
sx q[3];
rz(-1.2602826) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.96596232) q[0];
sx q[0];
rz(-1.9251134) q[0];
sx q[0];
rz(-2.8438582) q[0];
rz(-0.072602428) q[1];
sx q[1];
rz(-2.4563792) q[1];
sx q[1];
rz(0.6822449) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84083623) q[0];
sx q[0];
rz(-1.8634691) q[0];
sx q[0];
rz(-2.2038644) q[0];
rz(-1.7715681) q[2];
sx q[2];
rz(-0.64324311) q[2];
sx q[2];
rz(0.40293504) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.025921019) q[1];
sx q[1];
rz(-2.762315) q[1];
sx q[1];
rz(1.9801407) q[1];
x q[2];
rz(-0.37118427) q[3];
sx q[3];
rz(-2.1371578) q[3];
sx q[3];
rz(1.0584099) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.5630774) q[2];
sx q[2];
rz(-2.7882521) q[2];
sx q[2];
rz(-3.0541218) q[2];
rz(-0.89720094) q[3];
sx q[3];
rz(-1.8950491) q[3];
sx q[3];
rz(0.38952601) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0422269) q[0];
sx q[0];
rz(-2.0392188) q[0];
sx q[0];
rz(1.9788096) q[0];
rz(-0.21576628) q[1];
sx q[1];
rz(-2.3240604) q[1];
sx q[1];
rz(-2.20631) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1638086) q[0];
sx q[0];
rz(-0.80923128) q[0];
sx q[0];
rz(-2.8092395) q[0];
rz(-pi) q[1];
x q[1];
rz(2.9940656) q[2];
sx q[2];
rz(-1.8263683) q[2];
sx q[2];
rz(2.7856261) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7120119) q[1];
sx q[1];
rz(-1.4918527) q[1];
sx q[1];
rz(-0.78671766) q[1];
x q[2];
rz(1.2452081) q[3];
sx q[3];
rz(-1.2712443) q[3];
sx q[3];
rz(1.3249782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0995348) q[2];
sx q[2];
rz(-0.65379405) q[2];
sx q[2];
rz(0.81134861) q[2];
rz(-0.30284303) q[3];
sx q[3];
rz(-1.9767913) q[3];
sx q[3];
rz(0.6306878) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.52043668) q[0];
sx q[0];
rz(-0.32670894) q[0];
sx q[0];
rz(2.4503571) q[0];
rz(0.65903819) q[1];
sx q[1];
rz(-1.5182511) q[1];
sx q[1];
rz(-0.59115994) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49159971) q[0];
sx q[0];
rz(-0.2831471) q[0];
sx q[0];
rz(1.6121907) q[0];
x q[1];
rz(-2.8560258) q[2];
sx q[2];
rz(-2.0519678) q[2];
sx q[2];
rz(-0.27142957) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4860079) q[1];
sx q[1];
rz(-2.297245) q[1];
sx q[1];
rz(2.6996711) q[1];
x q[2];
rz(-1.3839327) q[3];
sx q[3];
rz(-2.8442743) q[3];
sx q[3];
rz(-0.0082810915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9968694) q[2];
sx q[2];
rz(-2.150712) q[2];
sx q[2];
rz(-1.4608176) q[2];
rz(-0.75374323) q[3];
sx q[3];
rz(-1.4494579) q[3];
sx q[3];
rz(0.47372216) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4562456) q[0];
sx q[0];
rz(-2.0401177) q[0];
sx q[0];
rz(-0.22612485) q[0];
rz(1.4489168) q[1];
sx q[1];
rz(-2.1782404) q[1];
sx q[1];
rz(-2.1400145) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6152214) q[0];
sx q[0];
rz(-0.60524056) q[0];
sx q[0];
rz(-2.8855118) q[0];
x q[1];
rz(1.1923157) q[2];
sx q[2];
rz(-0.90736249) q[2];
sx q[2];
rz(2.0326234) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.16538936) q[1];
sx q[1];
rz(-1.7586629) q[1];
sx q[1];
rz(-1.4632744) q[1];
rz(-pi) q[2];
rz(2.4572152) q[3];
sx q[3];
rz(-1.3357786) q[3];
sx q[3];
rz(-0.37684611) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.5292624) q[2];
sx q[2];
rz(-2.3713106) q[2];
sx q[2];
rz(-0.15753499) q[2];
rz(2.9869288) q[3];
sx q[3];
rz(-1.9779343) q[3];
sx q[3];
rz(-0.4304339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4100819) q[0];
sx q[0];
rz(-1.9022576) q[0];
sx q[0];
rz(-0.70247689) q[0];
rz(-0.53600535) q[1];
sx q[1];
rz(-0.82967007) q[1];
sx q[1];
rz(1.3943256) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6781021) q[0];
sx q[0];
rz(-0.16568389) q[0];
sx q[0];
rz(1.413762) q[0];
rz(-pi) q[1];
x q[1];
rz(2.2807806) q[2];
sx q[2];
rz(-2.99798) q[2];
sx q[2];
rz(-1.2705621) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.179217) q[1];
sx q[1];
rz(-2.0834647) q[1];
sx q[1];
rz(-2.8151399) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5709747) q[3];
sx q[3];
rz(-1.9863673) q[3];
sx q[3];
rz(-3.0316169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.006762) q[2];
sx q[2];
rz(-2.4595021) q[2];
sx q[2];
rz(1.4466059) q[2];
rz(0.26099482) q[3];
sx q[3];
rz(-2.1311396) q[3];
sx q[3];
rz(1.9058156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21988729) q[0];
sx q[0];
rz(-0.6548665) q[0];
sx q[0];
rz(-1.0153216) q[0];
rz(-2.0657397) q[1];
sx q[1];
rz(-1.8668108) q[1];
sx q[1];
rz(1.6783953) q[1];
rz(-1.4599316) q[2];
sx q[2];
rz(-1.9456429) q[2];
sx q[2];
rz(2.6425998) q[2];
rz(-0.76293972) q[3];
sx q[3];
rz(-1.8992103) q[3];
sx q[3];
rz(-2.7873399) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
