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
rz(-0.15260881) q[0];
sx q[0];
rz(-1.1996562) q[0];
sx q[0];
rz(-0.3056404) q[0];
rz(-2.8513554) q[1];
sx q[1];
rz(-1.4479535) q[1];
sx q[1];
rz(-1.0040959) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5619156) q[0];
sx q[0];
rz(-1.5439057) q[0];
sx q[0];
rz(-0.057232522) q[0];
rz(-1.8199241) q[2];
sx q[2];
rz(-1.3835356) q[2];
sx q[2];
rz(-1.934777) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.6184147) q[1];
sx q[1];
rz(-0.52408275) q[1];
sx q[1];
rz(0.89282234) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5283995) q[3];
sx q[3];
rz(-1.0323845) q[3];
sx q[3];
rz(1.9574894) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.64757887) q[2];
sx q[2];
rz(-0.61654377) q[2];
sx q[2];
rz(-0.54440633) q[2];
rz(-0.38862774) q[3];
sx q[3];
rz(-1.9184687) q[3];
sx q[3];
rz(-0.94038928) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.87043864) q[0];
sx q[0];
rz(-2.1156023) q[0];
sx q[0];
rz(1.1218659) q[0];
rz(2.0815966) q[1];
sx q[1];
rz(-2.6385939) q[1];
sx q[1];
rz(-0.88567919) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1141027) q[0];
sx q[0];
rz(-1.8654279) q[0];
sx q[0];
rz(0.31722169) q[0];
rz(-pi) q[1];
rz(0.76876872) q[2];
sx q[2];
rz(-0.93009206) q[2];
sx q[2];
rz(0.51962432) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.6180775) q[1];
sx q[1];
rz(-1.6517868) q[1];
sx q[1];
rz(-0.063132719) q[1];
rz(-pi) q[2];
rz(3.0360004) q[3];
sx q[3];
rz(-1.6882992) q[3];
sx q[3];
rz(1.8697704) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.9649967) q[2];
sx q[2];
rz(-1.4623888) q[2];
sx q[2];
rz(-3.1122567) q[2];
rz(2.7756179) q[3];
sx q[3];
rz(-0.47593203) q[3];
sx q[3];
rz(-2.9717145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40816864) q[0];
sx q[0];
rz(-1.3329196) q[0];
sx q[0];
rz(0.15342203) q[0];
rz(0.19639213) q[1];
sx q[1];
rz(-1.6729665) q[1];
sx q[1];
rz(-2.3207655) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3494075) q[0];
sx q[0];
rz(-2.2174412) q[0];
sx q[0];
rz(0.92766858) q[0];
rz(2.1961741) q[2];
sx q[2];
rz(-2.4074005) q[2];
sx q[2];
rz(1.9482833) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.094999639) q[1];
sx q[1];
rz(-0.41727704) q[1];
sx q[1];
rz(0.098578171) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7749805) q[3];
sx q[3];
rz(-1.5248581) q[3];
sx q[3];
rz(1.0339772) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.11744943) q[2];
sx q[2];
rz(-1.8154181) q[2];
sx q[2];
rz(0.93309012) q[2];
rz(-2.7945331) q[3];
sx q[3];
rz(-2.0324028) q[3];
sx q[3];
rz(0.6412653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1953122) q[0];
sx q[0];
rz(-1.1271789) q[0];
sx q[0];
rz(2.0830925) q[0];
rz(3.0463386) q[1];
sx q[1];
rz(-1.1493827) q[1];
sx q[1];
rz(-2.9375295) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6263555) q[0];
sx q[0];
rz(-1.7254819) q[0];
sx q[0];
rz(1.7912639) q[0];
rz(-pi) q[1];
rz(-2.0458895) q[2];
sx q[2];
rz(-1.1563325) q[2];
sx q[2];
rz(1.1865277) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.3151824) q[1];
sx q[1];
rz(-0.99821222) q[1];
sx q[1];
rz(1.3316283) q[1];
rz(-pi) q[2];
rz(1.7768562) q[3];
sx q[3];
rz(-2.2584174) q[3];
sx q[3];
rz(1.7182409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0846587) q[2];
sx q[2];
rz(-1.6997507) q[2];
sx q[2];
rz(1.320768) q[2];
rz(1.002958) q[3];
sx q[3];
rz(-1.2489677) q[3];
sx q[3];
rz(-0.5654208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4249307) q[0];
sx q[0];
rz(-1.2053763) q[0];
sx q[0];
rz(2.5004814) q[0];
rz(-0.59116108) q[1];
sx q[1];
rz(-1.4407651) q[1];
sx q[1];
rz(-1.444918) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4137694) q[0];
sx q[0];
rz(-1.1846847) q[0];
sx q[0];
rz(-1.5408433) q[0];
rz(-pi) q[1];
rz(1.314268) q[2];
sx q[2];
rz(-0.70665828) q[2];
sx q[2];
rz(-0.15308943) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.302204) q[1];
sx q[1];
rz(-2.6008718) q[1];
sx q[1];
rz(-2.6469321) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.0743477) q[3];
sx q[3];
rz(-0.58234057) q[3];
sx q[3];
rz(-2.989188) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.5833907) q[2];
sx q[2];
rz(-0.99339569) q[2];
sx q[2];
rz(0.85303419) q[2];
rz(2.9663626) q[3];
sx q[3];
rz(-1.8195189) q[3];
sx q[3];
rz(-0.62293735) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22195062) q[0];
sx q[0];
rz(-0.80606824) q[0];
sx q[0];
rz(-2.6958418) q[0];
rz(2.1469965) q[1];
sx q[1];
rz(-2.1371195) q[1];
sx q[1];
rz(1.4882784) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5827245) q[0];
sx q[0];
rz(-1.7877401) q[0];
sx q[0];
rz(-2.4198351) q[0];
rz(0.044537466) q[2];
sx q[2];
rz(-1.4967398) q[2];
sx q[2];
rz(-0.088852126) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.46785746) q[1];
sx q[1];
rz(-2.7731491) q[1];
sx q[1];
rz(2.175827) q[1];
rz(1.8976977) q[3];
sx q[3];
rz(-1.4964678) q[3];
sx q[3];
rz(1.2352218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.17066869) q[2];
sx q[2];
rz(-1.8985775) q[2];
sx q[2];
rz(2.1155913) q[2];
rz(-0.79801997) q[3];
sx q[3];
rz(-2.3012216) q[3];
sx q[3];
rz(-0.26871267) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2949424) q[0];
sx q[0];
rz(-1.7468528) q[0];
sx q[0];
rz(-2.7251439) q[0];
rz(1.7448447) q[1];
sx q[1];
rz(-1.5301219) q[1];
sx q[1];
rz(-1.6646615) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.653841) q[0];
sx q[0];
rz(-1.6076822) q[0];
sx q[0];
rz(-1.0168309) q[0];
rz(-pi) q[1];
rz(1.9878108) q[2];
sx q[2];
rz(-1.872135) q[2];
sx q[2];
rz(0.84433198) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1394964) q[1];
sx q[1];
rz(-2.0130035) q[1];
sx q[1];
rz(-0.79331974) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4879513) q[3];
sx q[3];
rz(-1.2330862) q[3];
sx q[3];
rz(-3.0520671) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.2485409) q[2];
sx q[2];
rz(-0.38022843) q[2];
sx q[2];
rz(-2.5377972) q[2];
rz(-2.2327312) q[3];
sx q[3];
rz(-1.4330319) q[3];
sx q[3];
rz(-2.042167) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.73868442) q[0];
sx q[0];
rz(-1.383902) q[0];
sx q[0];
rz(2.4201194) q[0];
rz(-1.8062704) q[1];
sx q[1];
rz(-2.0029009) q[1];
sx q[1];
rz(-1.2566176) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8523113) q[0];
sx q[0];
rz(-0.86332488) q[0];
sx q[0];
rz(-1.1723164) q[0];
rz(1.0851791) q[2];
sx q[2];
rz(-2.0943421) q[2];
sx q[2];
rz(-2.8741037) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.9049038) q[1];
sx q[1];
rz(-1.3373378) q[1];
sx q[1];
rz(1.1462565) q[1];
x q[2];
rz(-1.6001892) q[3];
sx q[3];
rz(-1.6163858) q[3];
sx q[3];
rz(0.33908333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1544372) q[2];
sx q[2];
rz(-0.64266959) q[2];
sx q[2];
rz(1.8670234) q[2];
rz(-2.4529723) q[3];
sx q[3];
rz(-1.6504811) q[3];
sx q[3];
rz(-0.0037746599) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0140728) q[0];
sx q[0];
rz(-1.2496244) q[0];
sx q[0];
rz(2.4697812) q[0];
rz(-1.5348684) q[1];
sx q[1];
rz(-2.913919) q[1];
sx q[1];
rz(2.6188376) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5005258) q[0];
sx q[0];
rz(-2.4198774) q[0];
sx q[0];
rz(-2.7428377) q[0];
x q[1];
rz(-1.9143064) q[2];
sx q[2];
rz(-1.6186423) q[2];
sx q[2];
rz(0.44105083) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.46314374) q[1];
sx q[1];
rz(-1.9903954) q[1];
sx q[1];
rz(1.1518008) q[1];
rz(-1.3626484) q[3];
sx q[3];
rz(-1.2537595) q[3];
sx q[3];
rz(-0.014118346) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8911924) q[2];
sx q[2];
rz(-2.4947391) q[2];
sx q[2];
rz(2.7817173) q[2];
rz(-1.322809) q[3];
sx q[3];
rz(-1.4578994) q[3];
sx q[3];
rz(-0.047164269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5322402) q[0];
sx q[0];
rz(-0.9541963) q[0];
sx q[0];
rz(0.94616079) q[0];
rz(3.1006475) q[1];
sx q[1];
rz(-1.2766726) q[1];
sx q[1];
rz(1.2650222) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8403137) q[0];
sx q[0];
rz(-1.095533) q[0];
sx q[0];
rz(-1.454157) q[0];
rz(-pi) q[1];
rz(-0.47750116) q[2];
sx q[2];
rz(-0.5029808) q[2];
sx q[2];
rz(1.8790085) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.62502669) q[1];
sx q[1];
rz(-1.0426765) q[1];
sx q[1];
rz(0.81104802) q[1];
rz(-pi) q[2];
x q[2];
rz(1.446883) q[3];
sx q[3];
rz(-0.37624761) q[3];
sx q[3];
rz(2.8543684) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1653183) q[2];
sx q[2];
rz(-1.9860622) q[2];
sx q[2];
rz(0.18132845) q[2];
rz(2.2881962) q[3];
sx q[3];
rz(-2.3386164) q[3];
sx q[3];
rz(1.3388504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4637909) q[0];
sx q[0];
rz(-1.3297357) q[0];
sx q[0];
rz(1.7672675) q[0];
rz(-1.455066) q[1];
sx q[1];
rz(-1.68119) q[1];
sx q[1];
rz(2.841058) q[1];
rz(-0.64114943) q[2];
sx q[2];
rz(-1.9998301) q[2];
sx q[2];
rz(2.395973) q[2];
rz(1.3131014) q[3];
sx q[3];
rz(-1.292376) q[3];
sx q[3];
rz(-1.8348909) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
