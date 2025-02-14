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
rz(2.1514819) q[0];
sx q[0];
rz(-0.70073849) q[0];
sx q[0];
rz(-0.75779688) q[0];
rz(-2.4611729) q[1];
sx q[1];
rz(-1.0935723) q[1];
sx q[1];
rz(3.0419066) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2282277) q[0];
sx q[0];
rz(-2.1036357) q[0];
sx q[0];
rz(-1.9060978) q[0];
rz(-2.6131479) q[2];
sx q[2];
rz(-0.64982729) q[2];
sx q[2];
rz(-0.53534269) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0079818) q[1];
sx q[1];
rz(-0.33002285) q[1];
sx q[1];
rz(-2.2225478) q[1];
rz(-pi) q[2];
rz(1.7277579) q[3];
sx q[3];
rz(-2.779134) q[3];
sx q[3];
rz(-0.80998224) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.0199355) q[2];
sx q[2];
rz(-1.1172373) q[2];
sx q[2];
rz(-0.026148671) q[2];
rz(1.7805685) q[3];
sx q[3];
rz(-1.0739505) q[3];
sx q[3];
rz(1.1092383) q[3];
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
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.015198135) q[0];
sx q[0];
rz(-2.4488738) q[0];
sx q[0];
rz(-1.7607081) q[0];
rz(-0.2805925) q[1];
sx q[1];
rz(-0.81604373) q[1];
sx q[1];
rz(2.7516344) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54444194) q[0];
sx q[0];
rz(-1.6100804) q[0];
sx q[0];
rz(0.76048373) q[0];
rz(3.105279) q[2];
sx q[2];
rz(-2.2507689) q[2];
sx q[2];
rz(0.84062409) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.38052032) q[1];
sx q[1];
rz(-1.7514844) q[1];
sx q[1];
rz(1.5086345) q[1];
x q[2];
rz(-2.6282464) q[3];
sx q[3];
rz(-1.5905375) q[3];
sx q[3];
rz(-1.4743559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.0754764) q[2];
sx q[2];
rz(-2.7710997) q[2];
sx q[2];
rz(-2.9541435) q[2];
rz(-2.1665393) q[3];
sx q[3];
rz(-1.2504028) q[3];
sx q[3];
rz(1.8072849) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27075574) q[0];
sx q[0];
rz(-0.90553951) q[0];
sx q[0];
rz(-1.4181197) q[0];
rz(2.1899147) q[1];
sx q[1];
rz(-0.27351174) q[1];
sx q[1];
rz(0.66114122) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1456297) q[0];
sx q[0];
rz(-2.0858156) q[0];
sx q[0];
rz(-1.4670232) q[0];
x q[1];
rz(0.44936835) q[2];
sx q[2];
rz(-0.74715464) q[2];
sx q[2];
rz(-2.4408057) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-0.97296827) q[1];
sx q[1];
rz(-1.428471) q[1];
sx q[1];
rz(-0.032802805) q[1];
x q[2];
rz(1.4641573) q[3];
sx q[3];
rz(-1.0141981) q[3];
sx q[3];
rz(2.9857991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.4989138) q[2];
sx q[2];
rz(-1.4475409) q[2];
sx q[2];
rz(-1.0399237) q[2];
rz(0.47636473) q[3];
sx q[3];
rz(-0.56291181) q[3];
sx q[3];
rz(-2.3001455) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.3570525) q[0];
sx q[0];
rz(-2.3226876) q[0];
sx q[0];
rz(-0.99405974) q[0];
rz(-1.1114063) q[1];
sx q[1];
rz(-1.2329654) q[1];
sx q[1];
rz(0.19817373) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1678038) q[0];
sx q[0];
rz(-2.7884169) q[0];
sx q[0];
rz(0.26539452) q[0];
rz(0.54951602) q[2];
sx q[2];
rz(-1.9870111) q[2];
sx q[2];
rz(0.93577318) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.930814) q[1];
sx q[1];
rz(-0.77180201) q[1];
sx q[1];
rz(-0.43369689) q[1];
rz(3.1311225) q[3];
sx q[3];
rz(-1.0426077) q[3];
sx q[3];
rz(1.2658324) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.17890945) q[2];
sx q[2];
rz(-1.6390025) q[2];
sx q[2];
rz(-0.34509125) q[2];
rz(1.0382477) q[3];
sx q[3];
rz(-1.1188666) q[3];
sx q[3];
rz(-0.071917608) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2800696) q[0];
sx q[0];
rz(-0.11174209) q[0];
sx q[0];
rz(0.34568632) q[0];
rz(3.0142036) q[1];
sx q[1];
rz(-2.7368339) q[1];
sx q[1];
rz(1.60188) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4152849) q[0];
sx q[0];
rz(-1.3506225) q[0];
sx q[0];
rz(2.9998695) q[0];
x q[1];
rz(0.56371828) q[2];
sx q[2];
rz(-2.3968389) q[2];
sx q[2];
rz(0.60840397) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.3469763) q[1];
sx q[1];
rz(-1.2656175) q[1];
sx q[1];
rz(0.48726122) q[1];
rz(-pi) q[2];
x q[2];
rz(0.13700142) q[3];
sx q[3];
rz(-2.5638201) q[3];
sx q[3];
rz(-0.96891257) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.6973528) q[2];
sx q[2];
rz(-1.3109861) q[2];
sx q[2];
rz(0.41651192) q[2];
rz(-3.0537649) q[3];
sx q[3];
rz(-2.646324) q[3];
sx q[3];
rz(2.7145568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.86513615) q[0];
sx q[0];
rz(-0.75053954) q[0];
sx q[0];
rz(0.99391341) q[0];
rz(1.3424501) q[1];
sx q[1];
rz(-1.3748704) q[1];
sx q[1];
rz(1.3329175) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48644984) q[0];
sx q[0];
rz(-0.37715366) q[0];
sx q[0];
rz(-2.7906453) q[0];
rz(-pi) q[1];
rz(1.3424043) q[2];
sx q[2];
rz(-1.9740385) q[2];
sx q[2];
rz(0.58951145) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.87403395) q[1];
sx q[1];
rz(-1.4324974) q[1];
sx q[1];
rz(-3.1048246) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.366863) q[3];
sx q[3];
rz(-1.9074554) q[3];
sx q[3];
rz(-0.62420732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.6028613) q[2];
sx q[2];
rz(-1.8659325) q[2];
sx q[2];
rz(-2.4219647) q[2];
rz(2.1384278) q[3];
sx q[3];
rz(-1.4041308) q[3];
sx q[3];
rz(-0.42705944) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7506943) q[0];
sx q[0];
rz(-2.7321132) q[0];
sx q[0];
rz(0.77792186) q[0];
rz(-1.3748417) q[1];
sx q[1];
rz(-0.96343103) q[1];
sx q[1];
rz(2.8533997) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8223234) q[0];
sx q[0];
rz(-2.6900969) q[0];
sx q[0];
rz(1.176693) q[0];
rz(0.32518816) q[2];
sx q[2];
rz(-3.1085494) q[2];
sx q[2];
rz(-0.059800241) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.3871875) q[1];
sx q[1];
rz(-1.8251694) q[1];
sx q[1];
rz(-0.53577975) q[1];
rz(-1.1918357) q[3];
sx q[3];
rz(-1.6917479) q[3];
sx q[3];
rz(2.9658776) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.9162743) q[2];
sx q[2];
rz(-0.83471966) q[2];
sx q[2];
rz(1.6342596) q[2];
rz(-0.39089125) q[3];
sx q[3];
rz(-1.8796128) q[3];
sx q[3];
rz(1.619537) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7759719) q[0];
sx q[0];
rz(-2.2052152) q[0];
sx q[0];
rz(0.17882624) q[0];
rz(1.7504182) q[1];
sx q[1];
rz(-1.4427253) q[1];
sx q[1];
rz(-0.19405445) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.80668215) q[0];
sx q[0];
rz(-1.1977949) q[0];
sx q[0];
rz(-1.0399517) q[0];
rz(-pi) q[1];
rz(-2.4670635) q[2];
sx q[2];
rz(-0.87376587) q[2];
sx q[2];
rz(-1.7473237) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0240101) q[1];
sx q[1];
rz(-0.76723209) q[1];
sx q[1];
rz(-1.3831286) q[1];
x q[2];
rz(-0.18004288) q[3];
sx q[3];
rz(-0.92623912) q[3];
sx q[3];
rz(2.8839437) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.27851823) q[2];
sx q[2];
rz(-0.14894177) q[2];
sx q[2];
rz(-1.8982559) q[2];
rz(-2.8362823) q[3];
sx q[3];
rz(-1.52933) q[3];
sx q[3];
rz(2.7436411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0447926) q[0];
sx q[0];
rz(-3.0718006) q[0];
sx q[0];
rz(3.0332562) q[0];
rz(0.70016742) q[1];
sx q[1];
rz(-1.533875) q[1];
sx q[1];
rz(0.38665032) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2359206) q[0];
sx q[0];
rz(-0.86109829) q[0];
sx q[0];
rz(-0.88881148) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.79603802) q[2];
sx q[2];
rz(-0.21547367) q[2];
sx q[2];
rz(-0.79677671) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.3503482) q[1];
sx q[1];
rz(-1.2859259) q[1];
sx q[1];
rz(0.71276748) q[1];
rz(-pi) q[2];
rz(0.4824494) q[3];
sx q[3];
rz(-1.0257991) q[3];
sx q[3];
rz(-1.6703796) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.75006968) q[2];
sx q[2];
rz(-1.1015026) q[2];
sx q[2];
rz(0.92292419) q[2];
rz(2.7643438) q[3];
sx q[3];
rz(-0.96830761) q[3];
sx q[3];
rz(1.7678461) q[3];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7467576) q[0];
sx q[0];
rz(-0.45811284) q[0];
sx q[0];
rz(-0.75570345) q[0];
rz(-2.34756) q[1];
sx q[1];
rz(-0.67818063) q[1];
sx q[1];
rz(-1.1620577) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.018022009) q[0];
sx q[0];
rz(-1.7777052) q[0];
sx q[0];
rz(1.1637972) q[0];
x q[1];
rz(1.1240578) q[2];
sx q[2];
rz(-0.96914712) q[2];
sx q[2];
rz(0.52858875) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.5291497) q[1];
sx q[1];
rz(-2.4172287) q[1];
sx q[1];
rz(-2.0606478) q[1];
rz(-pi) q[2];
rz(0.28323876) q[3];
sx q[3];
rz(-1.1822209) q[3];
sx q[3];
rz(0.76517867) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0843167) q[2];
sx q[2];
rz(-0.82547775) q[2];
sx q[2];
rz(-0.1262854) q[2];
rz(-0.59857357) q[3];
sx q[3];
rz(-1.6152265) q[3];
sx q[3];
rz(1.9954782) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79692688) q[0];
sx q[0];
rz(-1.666259) q[0];
sx q[0];
rz(1.689612) q[0];
rz(-1.7563734) q[1];
sx q[1];
rz(-1.9671556) q[1];
sx q[1];
rz(3.0671493) q[1];
rz(-0.66567265) q[2];
sx q[2];
rz(-0.72279983) q[2];
sx q[2];
rz(1.3574251) q[2];
rz(0.53530467) q[3];
sx q[3];
rz(-0.44620958) q[3];
sx q[3];
rz(-0.37887497) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
