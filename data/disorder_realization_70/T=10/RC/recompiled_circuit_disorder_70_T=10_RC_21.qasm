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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.91416042) q[0];
sx q[0];
rz(-2.3254546) q[0];
sx q[0];
rz(-0.039365191) q[0];
rz(-pi) q[1];
rz(-0.91648957) q[2];
sx q[2];
rz(-1.620703) q[2];
sx q[2];
rz(-2.8061342) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.0178535) q[1];
sx q[1];
rz(-0.98736073) q[1];
sx q[1];
rz(0.53637335) q[1];
rz(-pi) q[2];
rz(-3.0862023) q[3];
sx q[3];
rz(-2.7408544) q[3];
sx q[3];
rz(1.3415568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.9156076) q[2];
sx q[2];
rz(-0.51708022) q[2];
sx q[2];
rz(1.263164) q[2];
rz(-1.8566711) q[3];
sx q[3];
rz(-1.4611171) q[3];
sx q[3];
rz(-2.8485956) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9830575) q[0];
sx q[0];
rz(-1.8479269) q[0];
sx q[0];
rz(2.7040226) q[0];
rz(2.5105387) q[1];
sx q[1];
rz(-0.32866207) q[1];
sx q[1];
rz(2.9204869) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9843922) q[0];
sx q[0];
rz(-1.385014) q[0];
sx q[0];
rz(1.1093344) q[0];
x q[1];
rz(1.2731304) q[2];
sx q[2];
rz(-2.7781099) q[2];
sx q[2];
rz(-0.01089451) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.931103) q[1];
sx q[1];
rz(-2.0277129) q[1];
sx q[1];
rz(-1.4355852) q[1];
x q[2];
rz(0.17640555) q[3];
sx q[3];
rz(-2.6325912) q[3];
sx q[3];
rz(-0.41364663) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.21800403) q[2];
sx q[2];
rz(-1.4322586) q[2];
sx q[2];
rz(-0.58829266) q[2];
rz(-2.6925987) q[3];
sx q[3];
rz(-2.7136927) q[3];
sx q[3];
rz(2.0935521) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5730729) q[0];
sx q[0];
rz(-0.097671106) q[0];
sx q[0];
rz(-1.0004689) q[0];
rz(-0.72552848) q[1];
sx q[1];
rz(-1.0745445) q[1];
sx q[1];
rz(0.75769889) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.476186) q[0];
sx q[0];
rz(-1.2520257) q[0];
sx q[0];
rz(2.589059) q[0];
x q[1];
rz(-0.20747848) q[2];
sx q[2];
rz(-1.7277272) q[2];
sx q[2];
rz(0.66397882) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.1670907) q[1];
sx q[1];
rz(-2.3579512) q[1];
sx q[1];
rz(2.2553315) q[1];
x q[2];
rz(1.682231) q[3];
sx q[3];
rz(-2.206344) q[3];
sx q[3];
rz(-0.089126822) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.2086601) q[2];
sx q[2];
rz(-2.9186086) q[2];
sx q[2];
rz(-2.3051252) q[2];
rz(1.4423192) q[3];
sx q[3];
rz(-1.6675555) q[3];
sx q[3];
rz(-2.9377655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1214685) q[0];
sx q[0];
rz(-3.0614873) q[0];
sx q[0];
rz(0.51112038) q[0];
rz(3.064149) q[1];
sx q[1];
rz(-2.7192392) q[1];
sx q[1];
rz(1.3285332) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15645813) q[0];
sx q[0];
rz(-1.7256323) q[0];
sx q[0];
rz(-1.7865208) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6287597) q[2];
sx q[2];
rz(-0.83216681) q[2];
sx q[2];
rz(1.2733449) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.7718539) q[1];
sx q[1];
rz(-0.26680294) q[1];
sx q[1];
rz(1.5727444) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7926746) q[3];
sx q[3];
rz(-1.4639877) q[3];
sx q[3];
rz(-2.2282003) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.093734309) q[2];
sx q[2];
rz(-1.2183943) q[2];
sx q[2];
rz(-1.1070586) q[2];
rz(-0.47248653) q[3];
sx q[3];
rz(-1.5269591) q[3];
sx q[3];
rz(2.2842177) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1012786) q[0];
sx q[0];
rz(-2.2450228) q[0];
sx q[0];
rz(0.3381981) q[0];
rz(-1.8473373) q[1];
sx q[1];
rz(-1.5816403) q[1];
sx q[1];
rz(2.6370874) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0230334) q[0];
sx q[0];
rz(-2.1001864) q[0];
sx q[0];
rz(2.8708007) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1930824) q[2];
sx q[2];
rz(-1.1337122) q[2];
sx q[2];
rz(1.9405685) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.63197631) q[1];
sx q[1];
rz(-2.7512433) q[1];
sx q[1];
rz(-2.2401287) q[1];
rz(0.78318627) q[3];
sx q[3];
rz(-0.83055701) q[3];
sx q[3];
rz(0.98379788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.111104) q[2];
sx q[2];
rz(-0.8384853) q[2];
sx q[2];
rz(-1.6476691) q[2];
rz(1.4533639) q[3];
sx q[3];
rz(-2.2036392) q[3];
sx q[3];
rz(-1.3822314) q[3];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21022739) q[0];
sx q[0];
rz(-0.87451044) q[0];
sx q[0];
rz(2.0507623) q[0];
rz(-2.6018654) q[1];
sx q[1];
rz(-1.153774) q[1];
sx q[1];
rz(-2.9528023) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8180346) q[0];
sx q[0];
rz(-1.8206882) q[0];
sx q[0];
rz(-2.7701993) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2115057) q[2];
sx q[2];
rz(-2.852716) q[2];
sx q[2];
rz(1.8747683) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.4420538) q[1];
sx q[1];
rz(-0.25146723) q[1];
sx q[1];
rz(-0.26775189) q[1];
rz(-1.7536229) q[3];
sx q[3];
rz(-1.8086686) q[3];
sx q[3];
rz(-0.49926234) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.171689) q[2];
sx q[2];
rz(-2.1128555) q[2];
sx q[2];
rz(1.8264654) q[2];
rz(1.5054437) q[3];
sx q[3];
rz(-0.35741487) q[3];
sx q[3];
rz(1.2833387) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(-1.0514907) q[0];
sx q[0];
rz(-0.83054709) q[0];
sx q[0];
rz(-2.8175957) q[0];
rz(1.8404768) q[1];
sx q[1];
rz(-2.3062861) q[1];
sx q[1];
rz(-1.7623998) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1349072) q[0];
sx q[0];
rz(-2.8838257) q[0];
sx q[0];
rz(1.2447312) q[0];
rz(-pi) q[1];
rz(-2.6799913) q[2];
sx q[2];
rz(-1.201655) q[2];
sx q[2];
rz(2.387407) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0849689) q[1];
sx q[1];
rz(-1.8693722) q[1];
sx q[1];
rz(1.741239) q[1];
x q[2];
rz(0.84545387) q[3];
sx q[3];
rz(-2.3618556) q[3];
sx q[3];
rz(-0.23507915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.032701187) q[2];
sx q[2];
rz(-1.4543433) q[2];
sx q[2];
rz(-2.1984055) q[2];
rz(-0.33106783) q[3];
sx q[3];
rz(-1.3676684) q[3];
sx q[3];
rz(0.29512063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.11809764) q[0];
sx q[0];
rz(-0.58693111) q[0];
sx q[0];
rz(-2.3773637) q[0];
rz(-3.0006192) q[1];
sx q[1];
rz(-2.7455175) q[1];
sx q[1];
rz(1.8483298) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8338776) q[0];
sx q[0];
rz(-2.6617962) q[0];
sx q[0];
rz(-2.0108372) q[0];
rz(0.13883491) q[2];
sx q[2];
rz(-2.5359557) q[2];
sx q[2];
rz(-2.5790737) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5544719) q[1];
sx q[1];
rz(-1.9738102) q[1];
sx q[1];
rz(0.15092571) q[1];
rz(-pi) q[2];
x q[2];
rz(1.771365) q[3];
sx q[3];
rz(-1.1784394) q[3];
sx q[3];
rz(-1.4810824) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.3747037) q[2];
sx q[2];
rz(-1.039144) q[2];
sx q[2];
rz(-0.58132201) q[2];
rz(-0.86822048) q[3];
sx q[3];
rz(-0.41728443) q[3];
sx q[3];
rz(-1.3114312) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.54810369) q[0];
sx q[0];
rz(-0.5797525) q[0];
sx q[0];
rz(-1.8909489) q[0];
rz(2.1447694) q[1];
sx q[1];
rz(-2.2579028) q[1];
sx q[1];
rz(1.2876127) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.33289136) q[0];
sx q[0];
rz(-1.5960072) q[0];
sx q[0];
rz(-2.0602134) q[0];
rz(2.8224045) q[2];
sx q[2];
rz(-2.6337998) q[2];
sx q[2];
rz(0.14949456) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.77546706) q[1];
sx q[1];
rz(-0.93535103) q[1];
sx q[1];
rz(0.25823621) q[1];
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
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8210956) q[2];
sx q[2];
rz(-2.7536776) q[2];
sx q[2];
rz(-1.8851177) q[2];
rz(0.16658941) q[3];
sx q[3];
rz(-1.5766671) q[3];
sx q[3];
rz(-1.0725718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88084108) q[0];
sx q[0];
rz(-0.80417019) q[0];
sx q[0];
rz(2.8163731) q[0];
rz(1.1351769) q[1];
sx q[1];
rz(-2.4644641) q[1];
sx q[1];
rz(-2.7744055) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8242278) q[0];
sx q[0];
rz(-1.6019078) q[0];
sx q[0];
rz(3.0935862) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7634002) q[2];
sx q[2];
rz(-0.66236712) q[2];
sx q[2];
rz(-2.3364002) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.91017427) q[1];
sx q[1];
rz(-1.7703729) q[1];
sx q[1];
rz(0.66399666) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.48093421) q[3];
sx q[3];
rz(-2.1554865) q[3];
sx q[3];
rz(-1.1993053) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.24511589) q[2];
sx q[2];
rz(-2.3325236) q[2];
sx q[2];
rz(-2.0754576) q[2];
rz(-3.0623479) q[3];
sx q[3];
rz(-0.60278046) q[3];
sx q[3];
rz(-1.665303) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8463678) q[0];
sx q[0];
rz(-1.5548779) q[0];
sx q[0];
rz(1.5665733) q[0];
rz(2.8425343) q[1];
sx q[1];
rz(-0.60332861) q[1];
sx q[1];
rz(0.52287846) q[1];
rz(2.3921641) q[2];
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
