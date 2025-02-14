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
rz(-1.6406887) q[0];
sx q[0];
rz(-2.357643) q[0];
sx q[0];
rz(-3.0408903) q[0];
rz(1.7164224) q[1];
sx q[1];
rz(-2.4727688) q[1];
sx q[1];
rz(-1.5688904) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11450311) q[0];
sx q[0];
rz(-0.56124262) q[0];
sx q[0];
rz(2.9176913) q[0];
x q[1];
rz(-0.27218453) q[2];
sx q[2];
rz(-1.6707175) q[2];
sx q[2];
rz(-1.8109494) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.4266963) q[1];
sx q[1];
rz(-2.9116803) q[1];
sx q[1];
rz(-0.80674361) q[1];
rz(0.80683462) q[3];
sx q[3];
rz(-2.5111622) q[3];
sx q[3];
rz(-0.85496074) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0147741) q[2];
sx q[2];
rz(-2.7744881) q[2];
sx q[2];
rz(1.611562) q[2];
rz(-1.4989467) q[3];
sx q[3];
rz(-0.25224125) q[3];
sx q[3];
rz(2.4422395) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
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
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67795578) q[0];
sx q[0];
rz(-1.3709443) q[0];
sx q[0];
rz(-2.875705) q[0];
rz(-1.5222585) q[1];
sx q[1];
rz(-1.8306754) q[1];
sx q[1];
rz(-1.486385) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2280272) q[0];
sx q[0];
rz(-1.5782713) q[0];
sx q[0];
rz(0.77534239) q[0];
rz(-pi) q[1];
rz(-1.8114016) q[2];
sx q[2];
rz(-2.1212656) q[2];
sx q[2];
rz(-1.6110422) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(1.3656908) q[1];
sx q[1];
rz(-1.2925102) q[1];
sx q[1];
rz(-2.6551618) q[1];
rz(-1.099212) q[3];
sx q[3];
rz(-1.3986482) q[3];
sx q[3];
rz(-2.4936683) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.8365525) q[2];
sx q[2];
rz(-2.063648) q[2];
sx q[2];
rz(3.0744699) q[2];
rz(0.51131311) q[3];
sx q[3];
rz(-1.7917683) q[3];
sx q[3];
rz(2.2093723) q[3];
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
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2387282) q[0];
sx q[0];
rz(-1.6599864) q[0];
sx q[0];
rz(0.1486775) q[0];
rz(0.45064926) q[1];
sx q[1];
rz(-1.5573749) q[1];
sx q[1];
rz(-0.91948909) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7203569) q[0];
sx q[0];
rz(-1.3034058) q[0];
sx q[0];
rz(1.8084099) q[0];
rz(-pi) q[1];
x q[1];
rz(0.8834768) q[2];
sx q[2];
rz(-1.4191437) q[2];
sx q[2];
rz(-2.6014181) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.28314134) q[1];
sx q[1];
rz(-1.0325924) q[1];
sx q[1];
rz(1.0665757) q[1];
x q[2];
rz(0.50289233) q[3];
sx q[3];
rz(-0.38030312) q[3];
sx q[3];
rz(1.2740434) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.273169) q[2];
sx q[2];
rz(-2.0563075) q[2];
sx q[2];
rz(-2.2853509) q[2];
rz(-2.26561) q[3];
sx q[3];
rz(-2.7464505) q[3];
sx q[3];
rz(1.1538848) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3530537) q[0];
sx q[0];
rz(-1.5748698) q[0];
sx q[0];
rz(2.5153611) q[0];
rz(1.4405454) q[1];
sx q[1];
rz(-0.82653058) q[1];
sx q[1];
rz(2.4339719) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9191878) q[0];
sx q[0];
rz(-2.4116431) q[0];
sx q[0];
rz(1.6987503) q[0];
rz(-pi) q[1];
rz(-2.796299) q[2];
sx q[2];
rz(-2.3273988) q[2];
sx q[2];
rz(-1.0619927) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.6615078) q[1];
sx q[1];
rz(-2.8605177) q[1];
sx q[1];
rz(-2.8451189) q[1];
x q[2];
rz(-1.510101) q[3];
sx q[3];
rz(-2.6558721) q[3];
sx q[3];
rz(-0.39050366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.4746728) q[2];
sx q[2];
rz(-1.2771353) q[2];
sx q[2];
rz(-2.7478768) q[2];
rz(-2.3422824) q[3];
sx q[3];
rz(-2.7768713) q[3];
sx q[3];
rz(-1.5289615) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(1.2213152) q[0];
sx q[0];
rz(-0.83196139) q[0];
sx q[0];
rz(0.078201683) q[0];
rz(0.68866628) q[1];
sx q[1];
rz(-2.2024901) q[1];
sx q[1];
rz(2.2130373) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9720042) q[0];
sx q[0];
rz(-2.1336335) q[0];
sx q[0];
rz(1.9408976) q[0];
rz(-pi) q[1];
rz(1.9276728) q[2];
sx q[2];
rz(-1.8584826) q[2];
sx q[2];
rz(1.9299802) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.2714933) q[1];
sx q[1];
rz(-0.44259354) q[1];
sx q[1];
rz(2.6111994) q[1];
x q[2];
rz(-2.6070091) q[3];
sx q[3];
rz(-2.4494736) q[3];
sx q[3];
rz(2.2632556) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3277305) q[2];
sx q[2];
rz(-2.6321415) q[2];
sx q[2];
rz(-2.4895085) q[2];
rz(0.70513519) q[3];
sx q[3];
rz(-1.6465681) q[3];
sx q[3];
rz(0.42858538) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.240165) q[0];
sx q[0];
rz(-2.9225898) q[0];
sx q[0];
rz(1.3448311) q[0];
rz(-2.6118458) q[1];
sx q[1];
rz(-1.2836645) q[1];
sx q[1];
rz(1.3175499) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.12152762) q[0];
sx q[0];
rz(-2.1449267) q[0];
sx q[0];
rz(-0.0962538) q[0];
rz(-pi) q[1];
rz(1.3438247) q[2];
sx q[2];
rz(-1.3161873) q[2];
sx q[2];
rz(1.3019778) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.3287813) q[1];
sx q[1];
rz(-1.2320294) q[1];
sx q[1];
rz(0.12943204) q[1];
x q[2];
rz(-1.0617015) q[3];
sx q[3];
rz(-1.4056883) q[3];
sx q[3];
rz(0.77835876) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.162447) q[2];
sx q[2];
rz(-1.2787168) q[2];
sx q[2];
rz(-2.0549959) q[2];
rz(0.094001683) q[3];
sx q[3];
rz(-2.0407929) q[3];
sx q[3];
rz(0.47959685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.78733665) q[0];
sx q[0];
rz(-2.2945963) q[0];
sx q[0];
rz(0.081548668) q[0];
rz(-1.6185919) q[1];
sx q[1];
rz(-2.1269507) q[1];
sx q[1];
rz(0.60752216) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5089585) q[0];
sx q[0];
rz(-2.3130401) q[0];
sx q[0];
rz(-0.25018042) q[0];
rz(-pi) q[1];
rz(-0.46731575) q[2];
sx q[2];
rz(-1.9155587) q[2];
sx q[2];
rz(1.6925825) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.1178149) q[1];
sx q[1];
rz(-1.5573606) q[1];
sx q[1];
rz(1.4680844) q[1];
rz(-pi) q[2];
rz(-0.37346249) q[3];
sx q[3];
rz(-1.499553) q[3];
sx q[3];
rz(3.0687268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(3.0532694) q[2];
sx q[2];
rz(-2.93556) q[2];
sx q[2];
rz(0.44331178) q[2];
rz(0.045470227) q[3];
sx q[3];
rz(-1.9732405) q[3];
sx q[3];
rz(-0.62062353) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4745673) q[0];
sx q[0];
rz(-0.98282951) q[0];
sx q[0];
rz(-0.56812754) q[0];
rz(-1.8727632) q[1];
sx q[1];
rz(-1.1548235) q[1];
sx q[1];
rz(-2.5297129) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2308919) q[0];
sx q[0];
rz(-2.7534427) q[0];
sx q[0];
rz(-2.9064473) q[0];
rz(-pi) q[1];
rz(0.18576764) q[2];
sx q[2];
rz(-1.8427197) q[2];
sx q[2];
rz(-0.28803864) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.1398604) q[1];
sx q[1];
rz(-2.0563158) q[1];
sx q[1];
rz(-2.8667169) q[1];
rz(2.967784) q[3];
sx q[3];
rz(-1.7663071) q[3];
sx q[3];
rz(1.7988835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.2511217) q[2];
sx q[2];
rz(-1.8939563) q[2];
sx q[2];
rz(2.9337511) q[2];
rz(-0.017092997) q[3];
sx q[3];
rz(-1.0216252) q[3];
sx q[3];
rz(2.0489073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6297778) q[0];
sx q[0];
rz(-2.9371174) q[0];
sx q[0];
rz(-1.7275607) q[0];
rz(-3.0746225) q[1];
sx q[1];
rz(-1.8194865) q[1];
sx q[1];
rz(-3.0019143) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.84565879) q[0];
sx q[0];
rz(-2.686123) q[0];
sx q[0];
rz(-2.3607872) q[0];
rz(-pi) q[1];
x q[1];
rz(2.0099747) q[2];
sx q[2];
rz(-2.1358523) q[2];
sx q[2];
rz(-1.4840839) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5433257) q[1];
sx q[1];
rz(-1.7964592) q[1];
sx q[1];
rz(1.78701) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2788676) q[3];
sx q[3];
rz(-2.9649884) q[3];
sx q[3];
rz(-2.7677329) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.70860538) q[2];
sx q[2];
rz(-1.2689509) q[2];
sx q[2];
rz(-0.73680669) q[2];
rz(1.5271198) q[3];
sx q[3];
rz(-2.1325285) q[3];
sx q[3];
rz(-1.2705151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0438185) q[0];
sx q[0];
rz(-2.173824) q[0];
sx q[0];
rz(-1.9737825) q[0];
rz(0.66520005) q[1];
sx q[1];
rz(-1.8635112) q[1];
sx q[1];
rz(0.98242378) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5031122) q[0];
sx q[0];
rz(-1.3706722) q[0];
sx q[0];
rz(-1.0699231) q[0];
rz(-0.63913156) q[2];
sx q[2];
rz(-1.35372) q[2];
sx q[2];
rz(2.9242196) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.1878882) q[1];
sx q[1];
rz(-2.4361103) q[1];
sx q[1];
rz(-3.0445552) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1725964) q[3];
sx q[3];
rz(-1.6329582) q[3];
sx q[3];
rz(-2.8650995) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.3490225) q[2];
sx q[2];
rz(-2.4945365) q[2];
sx q[2];
rz(-1.6046074) q[2];
rz(1.0200621) q[3];
sx q[3];
rz(-1.6198817) q[3];
sx q[3];
rz(0.1608688) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1178602) q[0];
sx q[0];
rz(-1.5570138) q[0];
sx q[0];
rz(1.8099063) q[0];
rz(-0.82904077) q[1];
sx q[1];
rz(-1.6571028) q[1];
sx q[1];
rz(2.170457) q[1];
rz(2.7132158) q[2];
sx q[2];
rz(-0.72441341) q[2];
sx q[2];
rz(-2.7073467) q[2];
rz(1.0868514) q[3];
sx q[3];
rz(-1.1900569) q[3];
sx q[3];
rz(-2.2708683) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
