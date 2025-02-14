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
rz(-1.2054297) q[0];
sx q[0];
rz(-3.087145) q[0];
sx q[0];
rz(2.2757538) q[0];
rz(1.6021597) q[1];
sx q[1];
rz(-1.8973693) q[1];
sx q[1];
rz(3.08334) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4522432) q[0];
sx q[0];
rz(-1.6177243) q[0];
sx q[0];
rz(1.7353135) q[0];
rz(-pi) q[1];
x q[1];
rz(0.073926386) q[2];
sx q[2];
rz(-2.0652899) q[2];
sx q[2];
rz(1.489515) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.73535753) q[1];
sx q[1];
rz(-0.82338295) q[1];
sx q[1];
rz(2.3830692) q[1];
rz(-pi) q[2];
rz(-1.2236113) q[3];
sx q[3];
rz(-1.876374) q[3];
sx q[3];
rz(1.3936199) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.8286579) q[2];
sx q[2];
rz(-0.88064319) q[2];
sx q[2];
rz(0.9642967) q[2];
rz(-3.0621081) q[3];
sx q[3];
rz(-1.3026404) q[3];
sx q[3];
rz(-2.3704119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.013414772) q[0];
sx q[0];
rz(-2.7981813) q[0];
sx q[0];
rz(-1.5403904) q[0];
rz(0.84114289) q[1];
sx q[1];
rz(-1.4079739) q[1];
sx q[1];
rz(0.85743633) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4824366) q[0];
sx q[0];
rz(-1.1322339) q[0];
sx q[0];
rz(3.1355804) q[0];
x q[1];
rz(-0.964999) q[2];
sx q[2];
rz(-2.409796) q[2];
sx q[2];
rz(3.059405) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.7341566) q[1];
sx q[1];
rz(-0.70824558) q[1];
sx q[1];
rz(0.70869653) q[1];
rz(-pi) q[2];
rz(2.8174761) q[3];
sx q[3];
rz(-2.8480011) q[3];
sx q[3];
rz(-2.8568134) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.88605827) q[2];
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
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(2.7968314) q[0];
sx q[0];
rz(-1.826257) q[0];
sx q[0];
rz(-1.8094081) q[0];
rz(-1.4765129) q[1];
sx q[1];
rz(-1.3308728) q[1];
sx q[1];
rz(0.61980334) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7999511) q[0];
sx q[0];
rz(-1.5551621) q[0];
sx q[0];
rz(3.0642516) q[0];
rz(-pi) q[1];
rz(1.31525) q[2];
sx q[2];
rz(-1.7339306) q[2];
sx q[2];
rz(-2.852885) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7156652) q[1];
sx q[1];
rz(-1.2536548) q[1];
sx q[1];
rz(-1.6460544) q[1];
rz(-2.8768646) q[3];
sx q[3];
rz(-0.94784289) q[3];
sx q[3];
rz(2.5271551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-3.1406113) q[2];
sx q[2];
rz(-2.8162214) q[2];
sx q[2];
rz(-2.3532383) q[2];
rz(-1.2761448) q[3];
sx q[3];
rz(-1.2272464) q[3];
sx q[3];
rz(-2.2787794) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1514423) q[0];
sx q[0];
rz(-1.3862415) q[0];
sx q[0];
rz(-2.0748806) q[0];
rz(-1.7881296) q[1];
sx q[1];
rz(-2.2746494) q[1];
sx q[1];
rz(-1.1563168) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0321372) q[0];
sx q[0];
rz(-2.9988213) q[0];
sx q[0];
rz(-1.7032864) q[0];
rz(-pi) q[1];
rz(-0.99937352) q[2];
sx q[2];
rz(-0.048437645) q[2];
sx q[2];
rz(-2.9349309) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9741648) q[1];
sx q[1];
rz(-1.8596223) q[1];
sx q[1];
rz(-1.1812737) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.084454704) q[3];
sx q[3];
rz(-0.66612426) q[3];
sx q[3];
rz(-0.67557456) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.68916965) q[2];
sx q[2];
rz(-2.4781879) q[2];
sx q[2];
rz(-0.058825292) q[2];
rz(2.5022653) q[3];
sx q[3];
rz(-1.9202193) q[3];
sx q[3];
rz(-2.2842469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.095379742) q[0];
sx q[0];
rz(-1.7261427) q[0];
sx q[0];
rz(0.010490622) q[0];
rz(1.5785716) q[1];
sx q[1];
rz(-0.69067162) q[1];
sx q[1];
rz(0.48935997) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78016475) q[0];
sx q[0];
rz(-0.95236822) q[0];
sx q[0];
rz(1.7649107) q[0];
x q[1];
rz(2.0496558) q[2];
sx q[2];
rz(-0.74981028) q[2];
sx q[2];
rz(-2.9887226) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.3172878) q[1];
sx q[1];
rz(-2.2154059) q[1];
sx q[1];
rz(1.1238453) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5813835) q[3];
sx q[3];
rz(-1.6015617) q[3];
sx q[3];
rz(0.72584541) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(1.2995149) q[2];
sx q[2];
rz(-2.044951) q[2];
sx q[2];
rz(-0.00046029885) q[2];
rz(0.67356235) q[3];
sx q[3];
rz(-2.2431777) q[3];
sx q[3];
rz(2.4323997) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8507268) q[0];
sx q[0];
rz(-1.3141661) q[0];
sx q[0];
rz(1.8379743) q[0];
rz(-2.8438026) q[1];
sx q[1];
rz(-0.89556634) q[1];
sx q[1];
rz(2.8477125) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0336908) q[0];
sx q[0];
rz(-1.0009888) q[0];
sx q[0];
rz(-0.66500647) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.43115012) q[2];
sx q[2];
rz(-2.2047055) q[2];
sx q[2];
rz(-2.0946787) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.10295) q[1];
sx q[1];
rz(-1.0979196) q[1];
sx q[1];
rz(0.49943083) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.5081283) q[3];
sx q[3];
rz(-2.7567299) q[3];
sx q[3];
rz(-0.72609392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.6728354) q[2];
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
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4918168) q[0];
sx q[0];
rz(-0.62115541) q[0];
sx q[0];
rz(-2.5671) q[0];
rz(-0.54221398) q[1];
sx q[1];
rz(-1.503399) q[1];
sx q[1];
rz(0.39074674) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6764561) q[0];
sx q[0];
rz(-0.97873298) q[0];
sx q[0];
rz(1.3964064) q[0];
x q[1];
rz(1.4847267) q[2];
sx q[2];
rz(-2.59554) q[2];
sx q[2];
rz(2.4459185) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3022932) q[1];
sx q[1];
rz(-1.5652735) q[1];
sx q[1];
rz(-2.7248565) q[1];
rz(-1.5338906) q[3];
sx q[3];
rz(-1.9362984) q[3];
sx q[3];
rz(0.17084641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.6812402) q[2];
sx q[2];
rz(-1.1972903) q[2];
sx q[2];
rz(-0.61275068) q[2];
rz(0.98226205) q[3];
sx q[3];
rz(-1.3145612) q[3];
sx q[3];
rz(2.6764892) q[3];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9923582) q[0];
sx q[0];
rz(-1.5748011) q[0];
sx q[0];
rz(0.055334844) q[0];
rz(1.5570359) q[1];
sx q[1];
rz(-1.0880071) q[1];
sx q[1];
rz(-0.43992821) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2429139) q[0];
sx q[0];
rz(-1.1132332) q[0];
sx q[0];
rz(-1.8458864) q[0];
rz(-pi) q[1];
rz(2.6505846) q[2];
sx q[2];
rz(-1.7517118) q[2];
sx q[2];
rz(2.343585) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(3.098843) q[1];
sx q[1];
rz(-2.1196445) q[1];
sx q[1];
rz(1.9427981) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1371577) q[3];
sx q[3];
rz(-1.8105339) q[3];
sx q[3];
rz(2.9842003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.7679193) q[2];
sx q[2];
rz(-1.7934711) q[2];
sx q[2];
rz(-0.31201735) q[2];
rz(-0.9196552) q[3];
sx q[3];
rz(-1.7731526) q[3];
sx q[3];
rz(-2.9254204) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8298518) q[0];
sx q[0];
rz(-0.98167247) q[0];
sx q[0];
rz(3.1134636) q[0];
rz(-1.4460538) q[1];
sx q[1];
rz(-1.1808993) q[1];
sx q[1];
rz(0.45949724) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0579679) q[0];
sx q[0];
rz(-3.0959385) q[0];
sx q[0];
rz(2.245957) q[0];
x q[1];
rz(-2.3192603) q[2];
sx q[2];
rz(-1.8206545) q[2];
sx q[2];
rz(2.9248934) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.34950911) q[1];
sx q[1];
rz(-2.4448423) q[1];
sx q[1];
rz(1.3443483) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.13925456) q[3];
sx q[3];
rz(-1.1920658) q[3];
sx q[3];
rz(-1.9851353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.10099899) q[2];
sx q[2];
rz(-1.7969635) q[2];
sx q[2];
rz(2.8864268) q[2];
rz(-2.1549639) q[3];
sx q[3];
rz(-2.7268703) q[3];
sx q[3];
rz(1.423098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3808909) q[0];
sx q[0];
rz(-1.0699027) q[0];
sx q[0];
rz(-1.3421407) q[0];
rz(0.0019207151) q[1];
sx q[1];
rz(-2.0989959) q[1];
sx q[1];
rz(-1.049918) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7586138) q[0];
sx q[0];
rz(-0.26018026) q[0];
sx q[0];
rz(0.95917137) q[0];
rz(-1.7467612) q[2];
sx q[2];
rz(-1.5719885) q[2];
sx q[2];
rz(0.81874412) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.62419301) q[1];
sx q[1];
rz(-1.788793) q[1];
sx q[1];
rz(-2.7647074) q[1];
rz(1.330659) q[3];
sx q[3];
rz(-2.3970345) q[3];
sx q[3];
rz(1.9743686) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2125825) q[2];
sx q[2];
rz(-1.2139823) q[2];
sx q[2];
rz(-2.6455961) q[2];
rz(1.4604733) q[3];
sx q[3];
rz(-1.7998327) q[3];
sx q[3];
rz(1.9546485) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
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
rz(-2.7620508) q[0];
sx q[0];
rz(-2.0370146) q[0];
sx q[0];
rz(1.4203352) q[0];
rz(2.5024391) q[1];
sx q[1];
rz(-1.2624546) q[1];
sx q[1];
rz(0.75844567) q[1];
rz(-0.47472246) q[2];
sx q[2];
rz(-1.3480777) q[2];
sx q[2];
rz(2.011275) q[2];
rz(0.88981723) q[3];
sx q[3];
rz(-0.87796904) q[3];
sx q[3];
rz(1.8129391) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
