OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.45286173) q[0];
sx q[0];
rz(-0.17012574) q[0];
sx q[0];
rz(2.3556019) q[0];
rz(0.6056447) q[1];
sx q[1];
rz(-2.4906467) q[1];
sx q[1];
rz(-0.63408607) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2848628) q[0];
sx q[0];
rz(-2.3861109) q[0];
sx q[0];
rz(-1.5289686) q[0];
rz(0.91648957) q[2];
sx q[2];
rz(-1.620703) q[2];
sx q[2];
rz(2.8061342) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.69971985) q[1];
sx q[1];
rz(-0.77077121) q[1];
sx q[1];
rz(2.2295879) q[1];
x q[2];
rz(-3.0862023) q[3];
sx q[3];
rz(-2.7408544) q[3];
sx q[3];
rz(-1.8000359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2259851) q[2];
sx q[2];
rz(-0.51708022) q[2];
sx q[2];
rz(-1.8784286) q[2];
rz(-1.8566711) q[3];
sx q[3];
rz(-1.4611171) q[3];
sx q[3];
rz(0.29299709) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.15853515) q[0];
sx q[0];
rz(-1.8479269) q[0];
sx q[0];
rz(2.7040226) q[0];
rz(-2.5105387) q[1];
sx q[1];
rz(-2.8129306) q[1];
sx q[1];
rz(-0.22110573) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3723345) q[0];
sx q[0];
rz(-0.49494574) q[0];
sx q[0];
rz(-1.9702205) q[0];
rz(-pi) q[1];
x q[1];
rz(1.919585) q[2];
sx q[2];
rz(-1.4663327) q[2];
sx q[2];
rz(1.2806569) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-3.0529777) q[1];
sx q[1];
rz(-2.6664475) q[1];
sx q[1];
rz(-0.26762025) q[1];
x q[2];
rz(-0.50237327) q[3];
sx q[3];
rz(-1.6564192) q[3];
sx q[3];
rz(-1.3115713) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.9235886) q[2];
sx q[2];
rz(-1.709334) q[2];
sx q[2];
rz(2.5533) q[2];
rz(-0.44899392) q[3];
sx q[3];
rz(-2.7136927) q[3];
sx q[3];
rz(1.0480405) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56851971) q[0];
sx q[0];
rz(-0.097671106) q[0];
sx q[0];
rz(-2.1411238) q[0];
rz(2.4160642) q[1];
sx q[1];
rz(-1.0745445) q[1];
sx q[1];
rz(0.75769889) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0962778) q[0];
sx q[0];
rz(-1.0490388) q[0];
sx q[0];
rz(-1.9406712) q[0];
rz(-pi) q[1];
rz(2.4865815) q[2];
sx q[2];
rz(-0.25946028) q[2];
sx q[2];
rz(-2.8734145) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.0211027) q[1];
sx q[1];
rz(-1.1081401) q[1];
sx q[1];
rz(-2.2282269) q[1];
rz(-pi) q[2];
rz(-1.682231) q[3];
sx q[3];
rz(-2.206344) q[3];
sx q[3];
rz(-3.0524658) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.2086601) q[2];
sx q[2];
rz(-0.22298403) q[2];
sx q[2];
rz(0.83646742) q[2];
rz(-1.4423192) q[3];
sx q[3];
rz(-1.4740372) q[3];
sx q[3];
rz(-2.9377655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.02012415) q[0];
sx q[0];
rz(-0.080105372) q[0];
sx q[0];
rz(2.6304723) q[0];
rz(-0.077443667) q[1];
sx q[1];
rz(-0.42235342) q[1];
sx q[1];
rz(1.8130594) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15645813) q[0];
sx q[0];
rz(-1.4159604) q[0];
sx q[0];
rz(-1.3550718) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0765692) q[2];
sx q[2];
rz(-2.2708714) q[2];
sx q[2];
rz(1.9698524) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1991785) q[1];
sx q[1];
rz(-1.5713099) q[1];
sx q[1];
rz(1.8375988) q[1];
rz(-pi) q[2];
x q[2];
rz(3.0321211) q[3];
sx q[3];
rz(-1.7913892) q[3];
sx q[3];
rz(0.68144875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(3.0478583) q[2];
sx q[2];
rz(-1.2183943) q[2];
sx q[2];
rz(-2.034534) q[2];
rz(-2.6691061) q[3];
sx q[3];
rz(-1.6146336) q[3];
sx q[3];
rz(-0.85737491) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(3.1012786) q[0];
sx q[0];
rz(-0.89656985) q[0];
sx q[0];
rz(0.3381981) q[0];
rz(-1.2942554) q[1];
sx q[1];
rz(-1.5816403) q[1];
sx q[1];
rz(-2.6370874) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5915247) q[0];
sx q[0];
rz(-1.3378157) q[0];
sx q[0];
rz(2.1165119) q[0];
x q[1];
rz(-1.9485103) q[2];
sx q[2];
rz(-2.0078805) q[2];
sx q[2];
rz(1.2010241) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.5704755) q[1];
sx q[1];
rz(-1.332453) q[1];
sx q[1];
rz(-1.2586602) q[1];
rz(-pi) q[2];
x q[2];
rz(2.2284331) q[3];
sx q[3];
rz(-1.0201766) q[3];
sx q[3];
rz(-0.0084358128) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0304886) q[2];
sx q[2];
rz(-2.3031074) q[2];
sx q[2];
rz(-1.6476691) q[2];
rz(1.6882287) q[3];
sx q[3];
rz(-0.93795347) q[3];
sx q[3];
rz(1.7593613) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21022739) q[0];
sx q[0];
rz(-2.2670822) q[0];
sx q[0];
rz(-1.0908303) q[0];
rz(-0.53972721) q[1];
sx q[1];
rz(-1.9878186) q[1];
sx q[1];
rz(-2.9528023) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.990373) q[0];
sx q[0];
rz(-1.2114721) q[0];
sx q[0];
rz(-1.3034526) q[0];
rz(-pi) q[1];
rz(1.2115057) q[2];
sx q[2];
rz(-2.852716) q[2];
sx q[2];
rz(-1.2668244) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.69953883) q[1];
sx q[1];
rz(-2.8901254) q[1];
sx q[1];
rz(2.8738408) q[1];
rz(2.8998428) q[3];
sx q[3];
rz(-1.7484192) q[3];
sx q[3];
rz(2.1135981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.9699036) q[2];
sx q[2];
rz(-2.1128555) q[2];
sx q[2];
rz(1.3151273) q[2];
rz(1.6361489) q[3];
sx q[3];
rz(-2.7841778) q[3];
sx q[3];
rz(-1.858254) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.090102) q[0];
sx q[0];
rz(-0.83054709) q[0];
sx q[0];
rz(-2.8175957) q[0];
rz(1.3011159) q[1];
sx q[1];
rz(-0.83530656) q[1];
sx q[1];
rz(-1.7623998) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8934879) q[0];
sx q[0];
rz(-1.4890492) q[0];
sx q[0];
rz(1.8155314) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.71521476) q[2];
sx q[2];
rz(-0.5826125) q[2];
sx q[2];
rz(-1.4441393) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.0849689) q[1];
sx q[1];
rz(-1.2722204) q[1];
sx q[1];
rz(1.4003537) q[1];
rz(-pi) q[2];
rz(2.561065) q[3];
sx q[3];
rz(-2.1248098) q[3];
sx q[3];
rz(1.1298657) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(3.1088915) q[2];
sx q[2];
rz(-1.4543433) q[2];
sx q[2];
rz(-0.94318715) q[2];
rz(-2.8105248) q[3];
sx q[3];
rz(-1.3676684) q[3];
sx q[3];
rz(2.846472) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.023495) q[0];
sx q[0];
rz(-0.58693111) q[0];
sx q[0];
rz(-2.3773637) q[0];
rz(-0.14097342) q[1];
sx q[1];
rz(-0.39607513) q[1];
sx q[1];
rz(1.8483298) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4828669) q[0];
sx q[0];
rz(-1.7687161) q[0];
sx q[0];
rz(-2.0107962) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6663315) q[2];
sx q[2];
rz(-2.1697858) q[2];
sx q[2];
rz(-0.73087382) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.1845491) q[1];
sx q[1];
rz(-0.42889412) q[1];
sx q[1];
rz(-1.2317608) q[1];
rz(-pi) q[2];
rz(2.6929018) q[3];
sx q[3];
rz(-2.7033227) q[3];
sx q[3];
rz(1.9697619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3747037) q[2];
sx q[2];
rz(-1.039144) q[2];
sx q[2];
rz(0.58132201) q[2];
rz(0.86822048) q[3];
sx q[3];
rz(-2.7243082) q[3];
sx q[3];
rz(1.8301615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.54810369) q[0];
sx q[0];
rz(-0.5797525) q[0];
sx q[0];
rz(1.8909489) q[0];
rz(-0.99682322) q[1];
sx q[1];
rz(-0.88368982) q[1];
sx q[1];
rz(-1.2876127) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2851965) q[0];
sx q[0];
rz(-0.49001339) q[0];
sx q[0];
rz(1.6243837) q[0];
rz(-pi) q[1];
rz(-0.31918819) q[2];
sx q[2];
rz(-2.6337998) q[2];
sx q[2];
rz(0.14949456) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.3567644) q[1];
sx q[1];
rz(-0.67912662) q[1];
sx q[1];
rz(1.9041512) q[1];
rz(-pi) q[2];
x q[2];
rz(1.6279531) q[3];
sx q[3];
rz(-0.74712979) q[3];
sx q[3];
rz(-2.3627918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.320497) q[2];
sx q[2];
rz(-0.38791502) q[2];
sx q[2];
rz(1.8851177) q[2];
rz(0.16658941) q[3];
sx q[3];
rz(-1.5766671) q[3];
sx q[3];
rz(2.0690209) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2607516) q[0];
sx q[0];
rz(-2.3374225) q[0];
sx q[0];
rz(-0.32521954) q[0];
rz(2.0064158) q[1];
sx q[1];
rz(-2.4644641) q[1];
sx q[1];
rz(-0.36718711) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4627535) q[0];
sx q[0];
rz(-3.084393) q[0];
sx q[0];
rz(0.57533933) q[0];
rz(-2.514421) q[2];
sx q[2];
rz(-1.3417202) q[2];
sx q[2];
rz(2.0723745) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.2314184) q[1];
sx q[1];
rz(-1.7703729) q[1];
sx q[1];
rz(0.66399666) q[1];
x q[2];
rz(-2.1807947) q[3];
sx q[3];
rz(-2.4028117) q[3];
sx q[3];
rz(1.9566386) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8964768) q[2];
sx q[2];
rz(-2.3325236) q[2];
sx q[2];
rz(-1.0661351) q[2];
rz(-0.079244763) q[3];
sx q[3];
rz(-0.60278046) q[3];
sx q[3];
rz(-1.4762896) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29522482) q[0];
sx q[0];
rz(-1.5548779) q[0];
sx q[0];
rz(1.5665733) q[0];
rz(0.29905839) q[1];
sx q[1];
rz(-2.538264) q[1];
sx q[1];
rz(-2.6187142) q[1];
rz(-0.74942855) q[2];
sx q[2];
rz(-1.6288169) q[2];
sx q[2];
rz(-1.5734869) q[2];
rz(-1.3973665) q[3];
sx q[3];
rz(-2.4122824) q[3];
sx q[3];
rz(-1.1576049) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
