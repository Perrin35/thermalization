OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.8795348) q[0];
sx q[0];
rz(-1.4095925) q[0];
sx q[0];
rz(-1.4341266) q[0];
rz(0.6342451) q[1];
sx q[1];
rz(6.8847818) q[1];
sx q[1];
rz(9.8431982) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.7941064) q[0];
sx q[0];
rz(-0.54437629) q[0];
sx q[0];
rz(1.044859) q[0];
x q[1];
rz(-1.4346052) q[2];
sx q[2];
rz(-1.4601267) q[2];
sx q[2];
rz(-1.1510804) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.5601555) q[1];
sx q[1];
rz(-1.1168224) q[1];
sx q[1];
rz(-2.8698688) q[1];
x q[2];
rz(1.1316142) q[3];
sx q[3];
rz(-1.7712799) q[3];
sx q[3];
rz(2.2905614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.3216386) q[2];
sx q[2];
rz(-1.3691838) q[2];
sx q[2];
rz(-2.3036172) q[2];
rz(-2.6485802) q[3];
sx q[3];
rz(-0.27291441) q[3];
sx q[3];
rz(-3.0626007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-0.74719602) q[0];
sx q[0];
rz(-0.72421873) q[0];
sx q[0];
rz(1.2778506) q[0];
rz(-0.17678075) q[1];
sx q[1];
rz(-1.8272094) q[1];
sx q[1];
rz(2.7094254) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9958187) q[0];
sx q[0];
rz(-0.60129014) q[0];
sx q[0];
rz(-0.50253089) q[0];
x q[1];
rz(-0.400153) q[2];
sx q[2];
rz(-1.1740985) q[2];
sx q[2];
rz(0.19043365) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.6501112) q[1];
sx q[1];
rz(-1.1307798) q[1];
sx q[1];
rz(0.12932175) q[1];
rz(-3.1167332) q[3];
sx q[3];
rz(-2.1356574) q[3];
sx q[3];
rz(-0.75627518) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.54923487) q[2];
sx q[2];
rz(-1.2640415) q[2];
sx q[2];
rz(2.6548927) q[2];
rz(-1.3782079) q[3];
sx q[3];
rz(-1.2599726) q[3];
sx q[3];
rz(0.53282213) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.10087) q[0];
sx q[0];
rz(-2.3951055) q[0];
sx q[0];
rz(-0.41734636) q[0];
rz(1.6529282) q[1];
sx q[1];
rz(-0.54549837) q[1];
sx q[1];
rz(0.506385) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34844549) q[0];
sx q[0];
rz(-1.9217102) q[0];
sx q[0];
rz(-2.9532414) q[0];
rz(-pi) q[1];
rz(0.26489139) q[2];
sx q[2];
rz(-0.75220097) q[2];
sx q[2];
rz(-2.0843992) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.0326929) q[1];
sx q[1];
rz(-1.1316205) q[1];
sx q[1];
rz(1.9119309) q[1];
rz(1.2012464) q[3];
sx q[3];
rz(-1.7412211) q[3];
sx q[3];
rz(-0.91526645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.4425519) q[2];
sx q[2];
rz(-0.46135819) q[2];
sx q[2];
rz(2.55012) q[2];
rz(-2.5555723) q[3];
sx q[3];
rz(-1.2089217) q[3];
sx q[3];
rz(1.7104141) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9716924) q[0];
sx q[0];
rz(-2.5132892) q[0];
sx q[0];
rz(-2.3024094) q[0];
rz(3.1160141) q[1];
sx q[1];
rz(-2.4459116) q[1];
sx q[1];
rz(-1.5930088) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73491053) q[0];
sx q[0];
rz(-1.4797749) q[0];
sx q[0];
rz(-2.0216366) q[0];
x q[1];
rz(-0.84393878) q[2];
sx q[2];
rz(-0.91274777) q[2];
sx q[2];
rz(2.7968) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.81739391) q[1];
sx q[1];
rz(-2.3014268) q[1];
sx q[1];
rz(-2.6170931) q[1];
rz(-pi) q[2];
rz(-2.7110093) q[3];
sx q[3];
rz(-1.3123371) q[3];
sx q[3];
rz(0.57627288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8362391) q[2];
sx q[2];
rz(-0.88399115) q[2];
sx q[2];
rz(-3.0419066) q[2];
rz(-2.1827407) q[3];
sx q[3];
rz(-1.8226263) q[3];
sx q[3];
rz(1.3249741) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.56617671) q[0];
sx q[0];
rz(-1.382099) q[0];
sx q[0];
rz(2.8856522) q[0];
rz(-0.4610962) q[1];
sx q[1];
rz(-2.0979116) q[1];
sx q[1];
rz(2.3815313) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81239031) q[0];
sx q[0];
rz(-0.15604067) q[0];
sx q[0];
rz(-0.72326707) q[0];
rz(-pi) q[1];
rz(1.1561469) q[2];
sx q[2];
rz(-1.5125325) q[2];
sx q[2];
rz(2.0410048) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.6872245) q[1];
sx q[1];
rz(-1.4441274) q[1];
sx q[1];
rz(-1.6997937) q[1];
rz(2.2171668) q[3];
sx q[3];
rz(-1.2741718) q[3];
sx q[3];
rz(-0.38277205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.5710859) q[2];
sx q[2];
rz(-1.7692302) q[2];
sx q[2];
rz(2.467353) q[2];
rz(-2.9267866) q[3];
sx q[3];
rz(-0.45682296) q[3];
sx q[3];
rz(3.1242483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53133416) q[0];
sx q[0];
rz(-1.4704309) q[0];
sx q[0];
rz(1.9624788) q[0];
rz(-0.20482652) q[1];
sx q[1];
rz(-0.79524672) q[1];
sx q[1];
rz(-1.0669605) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5087591) q[0];
sx q[0];
rz(-1.5501223) q[0];
sx q[0];
rz(1.5563006) q[0];
rz(-pi) q[1];
rz(2.7166769) q[2];
sx q[2];
rz(-1.1669461) q[2];
sx q[2];
rz(2.0770819) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.31159196) q[1];
sx q[1];
rz(-0.82468669) q[1];
sx q[1];
rz(3.0768865) q[1];
rz(-2.3451869) q[3];
sx q[3];
rz(-1.4561597) q[3];
sx q[3];
rz(2.93626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.71010464) q[2];
sx q[2];
rz(-1.8523214) q[2];
sx q[2];
rz(-0.55994326) q[2];
rz(2.4152749) q[3];
sx q[3];
rz(-0.30877078) q[3];
sx q[3];
rz(2.8360951) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85754919) q[0];
sx q[0];
rz(-0.54656583) q[0];
sx q[0];
rz(1.42111) q[0];
rz(2.9395318) q[1];
sx q[1];
rz(-1.4338564) q[1];
sx q[1];
rz(0.85817671) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9439529) q[0];
sx q[0];
rz(-1.4240992) q[0];
sx q[0];
rz(-0.12978817) q[0];
x q[1];
rz(3.0965205) q[2];
sx q[2];
rz(-1.1606693) q[2];
sx q[2];
rz(-2.8571667) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.3837636) q[1];
sx q[1];
rz(-1.5997636) q[1];
sx q[1];
rz(-1.5860735) q[1];
rz(-pi) q[2];
rz(2.6303597) q[3];
sx q[3];
rz(-2.5066262) q[3];
sx q[3];
rz(-2.883203) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.8537366) q[2];
sx q[2];
rz(-2.6617472) q[2];
sx q[2];
rz(-1.8161592) q[2];
rz(-0.89007968) q[3];
sx q[3];
rz(-1.1471014) q[3];
sx q[3];
rz(2.1530698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6778075) q[0];
sx q[0];
rz(-2.5690434) q[0];
sx q[0];
rz(0.37471399) q[0];
rz(2.162714) q[1];
sx q[1];
rz(-2.4596877) q[1];
sx q[1];
rz(-1.3495061) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6357248) q[0];
sx q[0];
rz(-2.3946107) q[0];
sx q[0];
rz(-2.5542459) q[0];
rz(3.1153203) q[2];
sx q[2];
rz(-0.71825829) q[2];
sx q[2];
rz(2.9422613) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.6264682) q[1];
sx q[1];
rz(-1.5212458) q[1];
sx q[1];
rz(-1.7527761) q[1];
x q[2];
rz(1.3091062) q[3];
sx q[3];
rz(-1.040254) q[3];
sx q[3];
rz(1.4025276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4902041) q[2];
sx q[2];
rz(-0.75275246) q[2];
sx q[2];
rz(-2.6728969) q[2];
rz(-1.1941341) q[3];
sx q[3];
rz(-1.3041376) q[3];
sx q[3];
rz(-1.3635925) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
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
rz(-0.63012183) q[0];
sx q[0];
rz(-2.2352495) q[0];
sx q[0];
rz(1.2619031) q[0];
rz(-2.966554) q[1];
sx q[1];
rz(-1.1418399) q[1];
sx q[1];
rz(1.5375686) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81048548) q[0];
sx q[0];
rz(-1.2665505) q[0];
sx q[0];
rz(0.23181339) q[0];
rz(-1.6749009) q[2];
sx q[2];
rz(-1.9874007) q[2];
sx q[2];
rz(-2.5660851) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.21737145) q[1];
sx q[1];
rz(-1.4573759) q[1];
sx q[1];
rz(2.4356615) q[1];
rz(-1.1101301) q[3];
sx q[3];
rz(-1.2445645) q[3];
sx q[3];
rz(-2.5903451) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.039915446) q[2];
sx q[2];
rz(-2.1890409) q[2];
sx q[2];
rz(0.76134479) q[2];
rz(2.2411761) q[3];
sx q[3];
rz(-0.59949985) q[3];
sx q[3];
rz(-0.049023978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.802357) q[0];
sx q[0];
rz(-0.46827066) q[0];
sx q[0];
rz(-0.21690579) q[0];
rz(-0.63198173) q[1];
sx q[1];
rz(-1.4914373) q[1];
sx q[1];
rz(-0.95473081) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.023708658) q[0];
sx q[0];
rz(-2.3291991) q[0];
sx q[0];
rz(-2.6931767) q[0];
x q[1];
rz(-1.9913313) q[2];
sx q[2];
rz(-2.0805801) q[2];
sx q[2];
rz(-1.2812986) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.2659457) q[1];
sx q[1];
rz(-1.6469643) q[1];
sx q[1];
rz(-2.6850558) q[1];
x q[2];
rz(2.2050489) q[3];
sx q[3];
rz(-1.9359971) q[3];
sx q[3];
rz(-1.6983502) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.1251936) q[2];
sx q[2];
rz(-1.23896) q[2];
sx q[2];
rz(-0.94474244) q[2];
rz(-0.38481209) q[3];
sx q[3];
rz(-2.0231569) q[3];
sx q[3];
rz(0.95705664) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.64086296) q[0];
sx q[0];
rz(-0.50518112) q[0];
sx q[0];
rz(1.5541979) q[0];
rz(2.2568933) q[1];
sx q[1];
rz(-2.2365166) q[1];
sx q[1];
rz(2.8832163) q[1];
rz(2.2434071) q[2];
sx q[2];
rz(-0.85815103) q[2];
sx q[2];
rz(1.382538) q[2];
rz(1.1273884) q[3];
sx q[3];
rz(-1.756712) q[3];
sx q[3];
rz(-2.0851019) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
