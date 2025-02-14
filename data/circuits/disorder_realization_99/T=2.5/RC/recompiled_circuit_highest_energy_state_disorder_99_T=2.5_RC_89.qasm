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
rz(1.9361629) q[0];
sx q[0];
rz(3.087145) q[0];
sx q[0];
rz(8.5589391) q[0];
rz(1.6021597) q[1];
sx q[1];
rz(-1.8973693) q[1];
sx q[1];
rz(-0.058252637) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39393878) q[0];
sx q[0];
rz(-0.17102111) q[0];
sx q[0];
rz(-1.8500516) q[0];
rz(3.0676663) q[2];
sx q[2];
rz(-1.0763028) q[2];
sx q[2];
rz(1.489515) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.4076947) q[1];
sx q[1];
rz(-1.0419894) q[1];
sx q[1];
rz(-2.4771792) q[1];
rz(-2.8179161) q[3];
sx q[3];
rz(-1.9012682) q[3];
sx q[3];
rz(-3.0728473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.3129348) q[2];
sx q[2];
rz(-2.2609495) q[2];
sx q[2];
rz(-2.177296) q[2];
rz(3.0621081) q[3];
sx q[3];
rz(-1.8389523) q[3];
sx q[3];
rz(0.77118072) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.013414772) q[0];
sx q[0];
rz(-0.34341136) q[0];
sx q[0];
rz(1.6012023) q[0];
rz(-0.84114289) q[1];
sx q[1];
rz(-1.4079739) q[1];
sx q[1];
rz(-0.85743633) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6591561) q[0];
sx q[0];
rz(-1.1322339) q[0];
sx q[0];
rz(-3.1355804) q[0];
rz(-pi) q[1];
rz(-2.1765937) q[2];
sx q[2];
rz(-0.73179663) q[2];
sx q[2];
rz(-0.082187637) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.5799478) q[1];
sx q[1];
rz(-2.0873293) q[1];
sx q[1];
rz(-1.0622611) q[1];
rz(1.6667834) q[3];
sx q[3];
rz(-1.848683) q[3];
sx q[3];
rz(3.0888219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.2555344) q[2];
sx q[2];
rz(-0.24571358) q[2];
sx q[2];
rz(-1.9319755) q[2];
rz(1.7602734) q[3];
sx q[3];
rz(-1.8807024) q[3];
sx q[3];
rz(2.5513726) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34476122) q[0];
sx q[0];
rz(-1.826257) q[0];
sx q[0];
rz(1.3321846) q[0];
rz(-1.4765129) q[1];
sx q[1];
rz(-1.8107199) q[1];
sx q[1];
rz(-0.61980334) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7999511) q[0];
sx q[0];
rz(-1.5551621) q[0];
sx q[0];
rz(-3.0642516) q[0];
rz(-pi) q[1];
rz(-1.31525) q[2];
sx q[2];
rz(-1.7339306) q[2];
sx q[2];
rz(-0.28870764) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.42592749) q[1];
sx q[1];
rz(-1.2536548) q[1];
sx q[1];
rz(-1.6460544) q[1];
rz(2.8768646) q[3];
sx q[3];
rz(-2.1937498) q[3];
sx q[3];
rz(2.5271551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.0009813112) q[2];
sx q[2];
rz(-0.32537127) q[2];
sx q[2];
rz(0.78835431) q[2];
rz(1.8654478) q[3];
sx q[3];
rz(-1.9143462) q[3];
sx q[3];
rz(2.2787794) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1514423) q[0];
sx q[0];
rz(-1.7553512) q[0];
sx q[0];
rz(1.066712) q[0];
rz(-1.7881296) q[1];
sx q[1];
rz(-0.86694327) q[1];
sx q[1];
rz(1.1563168) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3301834) q[0];
sx q[0];
rz(-1.5519987) q[0];
sx q[0];
rz(1.4292595) q[0];
x q[1];
rz(-1.5300445) q[2];
sx q[2];
rz(-1.5446071) q[2];
sx q[2];
rz(-1.9350236) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.51976065) q[1];
sx q[1];
rz(-1.9433738) q[1];
sx q[1];
rz(-0.31078679) q[1];
x q[2];
rz(2.4772024) q[3];
sx q[3];
rz(-1.6229462) q[3];
sx q[3];
rz(-0.96168226) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.68916965) q[2];
sx q[2];
rz(-2.4781879) q[2];
sx q[2];
rz(-3.0827674) q[2];
rz(0.63932738) q[3];
sx q[3];
rz(-1.2213734) q[3];
sx q[3];
rz(-2.2842469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0462129) q[0];
sx q[0];
rz(-1.4154499) q[0];
sx q[0];
rz(3.131102) q[0];
rz(-1.5785716) q[1];
sx q[1];
rz(-2.450921) q[1];
sx q[1];
rz(0.48935997) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45324126) q[0];
sx q[0];
rz(-2.4972389) q[0];
sx q[0];
rz(-0.2647361) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2614865) q[2];
sx q[2];
rz(-1.2513796) q[2];
sx q[2];
rz(1.781014) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.64397821) q[1];
sx q[1];
rz(-0.76581875) q[1];
sx q[1];
rz(0.52180565) q[1];
x q[2];
rz(-0.33133026) q[3];
sx q[3];
rz(-0.032535527) q[3];
sx q[3];
rz(2.0842541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8420777) q[2];
sx q[2];
rz(-1.0966417) q[2];
sx q[2];
rz(0.00046029885) q[2];
rz(2.4680303) q[3];
sx q[3];
rz(-2.2431777) q[3];
sx q[3];
rz(0.7091929) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8507268) q[0];
sx q[0];
rz(-1.3141661) q[0];
sx q[0];
rz(1.3036183) q[0];
rz(0.29779008) q[1];
sx q[1];
rz(-2.2460263) q[1];
sx q[1];
rz(0.29388014) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0760982) q[0];
sx q[0];
rz(-0.84653234) q[0];
sx q[0];
rz(2.3374071) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.05386) q[2];
sx q[2];
rz(-2.392025) q[2];
sx q[2];
rz(0.38656879) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.1634632) q[1];
sx q[1];
rz(-0.67369622) q[1];
sx q[1];
rz(-0.81836318) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.31556337) q[3];
sx q[3];
rz(-1.3466918) q[3];
sx q[3];
rz(1.4423646) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.46875724) q[2];
sx q[2];
rz(-1.7054649) q[2];
sx q[2];
rz(0.22077665) q[2];
rz(1.8686434) q[3];
sx q[3];
rz(-1.7468529) q[3];
sx q[3];
rz(-1.9934191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4918168) q[0];
sx q[0];
rz(-2.5204372) q[0];
sx q[0];
rz(2.5671) q[0];
rz(-2.5993787) q[1];
sx q[1];
rz(-1.6381936) q[1];
sx q[1];
rz(0.39074674) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46513656) q[0];
sx q[0];
rz(-0.97873298) q[0];
sx q[0];
rz(1.7451863) q[0];
rz(-pi) q[1];
rz(-1.656866) q[2];
sx q[2];
rz(-0.54605267) q[2];
sx q[2];
rz(-2.4459185) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.83929944) q[1];
sx q[1];
rz(-1.5763192) q[1];
sx q[1];
rz(-0.41673613) q[1];
rz(2.7758632) q[3];
sx q[3];
rz(-1.6052632) q[3];
sx q[3];
rz(-1.3867539) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.4603525) q[2];
sx q[2];
rz(-1.1972903) q[2];
sx q[2];
rz(0.61275068) q[2];
rz(0.98226205) q[3];
sx q[3];
rz(-1.8270315) q[3];
sx q[3];
rz(-2.6764892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9923582) q[0];
sx q[0];
rz(-1.5667916) q[0];
sx q[0];
rz(0.055334844) q[0];
rz(-1.5570359) q[1];
sx q[1];
rz(-1.0880071) q[1];
sx q[1];
rz(-2.7016644) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9377559) q[0];
sx q[0];
rz(-1.3246264) q[0];
sx q[0];
rz(2.6686644) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7753168) q[2];
sx q[2];
rz(-1.0884945) q[2];
sx q[2];
rz(0.67687809) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.042749636) q[1];
sx q[1];
rz(-2.1196445) q[1];
sx q[1];
rz(1.9427981) q[1];
rz(-pi) q[2];
rz(-0.28195076) q[3];
sx q[3];
rz(-1.022517) q[3];
sx q[3];
rz(-1.5632526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3736734) q[2];
sx q[2];
rz(-1.3481216) q[2];
sx q[2];
rz(-2.8295753) q[2];
rz(-2.2219374) q[3];
sx q[3];
rz(-1.3684401) q[3];
sx q[3];
rz(-2.9254204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8298518) q[0];
sx q[0];
rz(-2.1599202) q[0];
sx q[0];
rz(-0.028129015) q[0];
rz(1.6955388) q[1];
sx q[1];
rz(-1.9606934) q[1];
sx q[1];
rz(2.6820954) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0579679) q[0];
sx q[0];
rz(-0.045654181) q[0];
sx q[0];
rz(2.245957) q[0];
x q[1];
rz(2.3192603) q[2];
sx q[2];
rz(-1.3209381) q[2];
sx q[2];
rz(2.9248934) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.64134899) q[1];
sx q[1];
rz(-0.8952221) q[1];
sx q[1];
rz(-2.9558923) q[1];
x q[2];
rz(-1.2351843) q[3];
sx q[3];
rz(-0.40235717) q[3];
sx q[3];
rz(1.5188007) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.10099899) q[2];
sx q[2];
rz(-1.3446292) q[2];
sx q[2];
rz(0.25516587) q[2];
rz(0.98662871) q[3];
sx q[3];
rz(-2.7268703) q[3];
sx q[3];
rz(-1.7184947) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7607018) q[0];
sx q[0];
rz(-2.07169) q[0];
sx q[0];
rz(-1.3421407) q[0];
rz(-3.1396719) q[1];
sx q[1];
rz(-2.0989959) q[1];
sx q[1];
rz(2.0916746) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5494019) q[0];
sx q[0];
rz(-1.4225385) q[0];
sx q[0];
rz(1.3561983) q[0];
x q[1];
rz(1.5776063) q[2];
sx q[2];
rz(-0.17596888) q[2];
sx q[2];
rz(2.3962452) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.109595) q[1];
sx q[1];
rz(-1.9383311) q[1];
sx q[1];
rz(1.8046741) q[1];
rz(-0.2157507) q[3];
sx q[3];
rz(-0.85236824) q[3];
sx q[3];
rz(1.6529447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.2125825) q[2];
sx q[2];
rz(-1.9276103) q[2];
sx q[2];
rz(-0.49599656) q[2];
rz(1.4604733) q[3];
sx q[3];
rz(-1.3417599) q[3];
sx q[3];
rz(1.1869441) q[3];
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
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7620508) q[0];
sx q[0];
rz(-2.0370146) q[0];
sx q[0];
rz(1.4203352) q[0];
rz(0.63915359) q[1];
sx q[1];
rz(-1.879138) q[1];
sx q[1];
rz(-2.383147) q[1];
rz(1.3214618) q[2];
sx q[2];
rz(-2.032866) q[2];
sx q[2];
rz(0.55351071) q[2];
rz(-0.88981723) q[3];
sx q[3];
rz(-2.2636236) q[3];
sx q[3];
rz(-1.3286535) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
