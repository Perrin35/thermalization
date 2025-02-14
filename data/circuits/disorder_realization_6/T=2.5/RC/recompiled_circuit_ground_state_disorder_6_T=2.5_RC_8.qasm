OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.71894574) q[0];
sx q[0];
rz(-2.5526241) q[0];
sx q[0];
rz(-0.11864057) q[0];
rz(-0.9912107) q[1];
sx q[1];
rz(-2.0999496) q[1];
sx q[1];
rz(-0.51377327) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2048747) q[0];
sx q[0];
rz(-0.74130171) q[0];
sx q[0];
rz(0.78420774) q[0];
rz(-pi) q[1];
x q[1];
rz(0.47494048) q[2];
sx q[2];
rz(-1.5154334) q[2];
sx q[2];
rz(-2.9621015) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.39756718) q[1];
sx q[1];
rz(-1.356522) q[1];
sx q[1];
rz(-0.32886966) q[1];
rz(-pi) q[2];
x q[2];
rz(0.051187201) q[3];
sx q[3];
rz(-2.5489106) q[3];
sx q[3];
rz(0.41646233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.4673246) q[2];
sx q[2];
rz(-1.6215723) q[2];
sx q[2];
rz(2.8094214) q[2];
rz(-0.25003555) q[3];
sx q[3];
rz(-1.2269521) q[3];
sx q[3];
rz(1.8488319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8914723) q[0];
sx q[0];
rz(-1.8880867) q[0];
sx q[0];
rz(2.6943595) q[0];
rz(-1.0294634) q[1];
sx q[1];
rz(-0.58164683) q[1];
sx q[1];
rz(-0.10118016) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3295022) q[0];
sx q[0];
rz(-2.0116099) q[0];
sx q[0];
rz(-2.0157218) q[0];
rz(1.871505) q[2];
sx q[2];
rz(-0.60023396) q[2];
sx q[2];
rz(-2.9832961) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.356368) q[1];
sx q[1];
rz(-1.1572946) q[1];
sx q[1];
rz(-1.861839) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4199791) q[3];
sx q[3];
rz(-2.4972417) q[3];
sx q[3];
rz(-0.86101139) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0826565) q[2];
sx q[2];
rz(-0.75581789) q[2];
sx q[2];
rz(1.6015046) q[2];
rz(2.5236409) q[3];
sx q[3];
rz(-0.58303419) q[3];
sx q[3];
rz(1.1624973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5471632) q[0];
sx q[0];
rz(-1.0048486) q[0];
sx q[0];
rz(2.6255703) q[0];
rz(-1.0524606) q[1];
sx q[1];
rz(-2.2155589) q[1];
sx q[1];
rz(-1.9693536) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7965423) q[0];
sx q[0];
rz(-1.5597412) q[0];
sx q[0];
rz(2.7175236) q[0];
rz(2.2812821) q[2];
sx q[2];
rz(-1.6261618) q[2];
sx q[2];
rz(-2.4982128) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-1.7683923) q[1];
sx q[1];
rz(-0.52881634) q[1];
sx q[1];
rz(-1.9463198) q[1];
rz(-1.5564136) q[3];
sx q[3];
rz(-0.77061117) q[3];
sx q[3];
rz(1.5387467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.8664794) q[2];
sx q[2];
rz(-1.2629513) q[2];
sx q[2];
rz(-0.67683721) q[2];
rz(2.6721241) q[3];
sx q[3];
rz(-2.2468086) q[3];
sx q[3];
rz(-1.2137871) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4645828) q[0];
sx q[0];
rz(-2.1362169) q[0];
sx q[0];
rz(1.9418035) q[0];
rz(-2.5489573) q[1];
sx q[1];
rz(-1.0721782) q[1];
sx q[1];
rz(-1.4592272) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1593247) q[0];
sx q[0];
rz(-1.766681) q[0];
sx q[0];
rz(-2.5513493) q[0];
x q[1];
rz(-2.3427354) q[2];
sx q[2];
rz(-1.7513382) q[2];
sx q[2];
rz(1.7910166) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.693599) q[1];
sx q[1];
rz(-1.3136615) q[1];
sx q[1];
rz(0.93235459) q[1];
x q[2];
rz(1.9497037) q[3];
sx q[3];
rz(-0.92273308) q[3];
sx q[3];
rz(-2.5398382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.6478641) q[2];
sx q[2];
rz(-1.0961327) q[2];
sx q[2];
rz(-2.4596821) q[2];
rz(-2.72825) q[3];
sx q[3];
rz(-1.6091434) q[3];
sx q[3];
rz(1.1558862) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2979564) q[0];
sx q[0];
rz(-1.6652668) q[0];
sx q[0];
rz(-0.27808878) q[0];
rz(-0.9043215) q[1];
sx q[1];
rz(-1.6878637) q[1];
sx q[1];
rz(-0.98731891) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.322418) q[0];
sx q[0];
rz(-1.8083982) q[0];
sx q[0];
rz(0.95652076) q[0];
rz(-pi) q[1];
rz(-1.5426719) q[2];
sx q[2];
rz(-1.4338981) q[2];
sx q[2];
rz(1.2281451) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(3.1121171) q[1];
sx q[1];
rz(-1.6320845) q[1];
sx q[1];
rz(-1.960683) q[1];
rz(-pi) q[2];
rz(-1.0342399) q[3];
sx q[3];
rz(-1.3349541) q[3];
sx q[3];
rz(-0.72808108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.014293369) q[2];
sx q[2];
rz(-0.96724808) q[2];
sx q[2];
rz(0.75931749) q[2];
rz(1.6589818) q[3];
sx q[3];
rz(-1.9448152) q[3];
sx q[3];
rz(-0.35879859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.21869126) q[0];
sx q[0];
rz(-0.8194812) q[0];
sx q[0];
rz(-2.4500093) q[0];
rz(-2.2870731) q[1];
sx q[1];
rz(-2.4311192) q[1];
sx q[1];
rz(-2.7202594) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2153361) q[0];
sx q[0];
rz(-0.43784446) q[0];
sx q[0];
rz(1.7329751) q[0];
rz(-pi) q[1];
rz(2.8745804) q[2];
sx q[2];
rz(-1.8521554) q[2];
sx q[2];
rz(-0.54131484) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-3.1189733) q[1];
sx q[1];
rz(-1.2105064) q[1];
sx q[1];
rz(1.6605571) q[1];
x q[2];
rz(-0.38479067) q[3];
sx q[3];
rz(-2.5005385) q[3];
sx q[3];
rz(-1.3271063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.2104346) q[2];
sx q[2];
rz(-1.4782108) q[2];
sx q[2];
rz(-0.20273905) q[2];
rz(1.5984795) q[3];
sx q[3];
rz(-0.32035443) q[3];
sx q[3];
rz(1.6725484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8572674) q[0];
sx q[0];
rz(-1.4284644) q[0];
sx q[0];
rz(-0.39886928) q[0];
rz(1.9786037) q[1];
sx q[1];
rz(-0.82610026) q[1];
sx q[1];
rz(-2.861048) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4568854) q[0];
sx q[0];
rz(-2.5734757) q[0];
sx q[0];
rz(2.3756308) q[0];
rz(-pi) q[1];
rz(2.0899827) q[2];
sx q[2];
rz(-2.3943442) q[2];
sx q[2];
rz(2.0488536) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.032256) q[1];
sx q[1];
rz(-1.2660626) q[1];
sx q[1];
rz(-2.017657) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1637502) q[3];
sx q[3];
rz(-2.2740318) q[3];
sx q[3];
rz(3.0271526) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3205388) q[2];
sx q[2];
rz(-2.1619449) q[2];
sx q[2];
rz(1.6542356) q[2];
rz(-2.1982543) q[3];
sx q[3];
rz(-0.92883674) q[3];
sx q[3];
rz(1.1613891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5281552) q[0];
sx q[0];
rz(-1.4035839) q[0];
sx q[0];
rz(-2.6392537) q[0];
rz(-2.595064) q[1];
sx q[1];
rz(-1.2556475) q[1];
sx q[1];
rz(-0.46195236) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97945178) q[0];
sx q[0];
rz(-2.5573423) q[0];
sx q[0];
rz(1.0229848) q[0];
rz(-pi) q[1];
rz(2.2221379) q[2];
sx q[2];
rz(-1.6088789) q[2];
sx q[2];
rz(3.0236261) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.5344055) q[1];
sx q[1];
rz(-2.7173373) q[1];
sx q[1];
rz(-1.8380182) q[1];
rz(-1.5835207) q[3];
sx q[3];
rz(-2.1012348) q[3];
sx q[3];
rz(0.96543559) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.168557) q[2];
sx q[2];
rz(-2.904197) q[2];
sx q[2];
rz(-0.17042223) q[2];
rz(-0.45525822) q[3];
sx q[3];
rz(-1.596343) q[3];
sx q[3];
rz(-2.4236603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1667204) q[0];
sx q[0];
rz(-1.6396739) q[0];
sx q[0];
rz(2.9611452) q[0];
rz(-0.51668733) q[1];
sx q[1];
rz(-1.5645212) q[1];
sx q[1];
rz(0.43389854) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1922197) q[0];
sx q[0];
rz(-1.9967955) q[0];
sx q[0];
rz(-2.1550687) q[0];
x q[1];
rz(-2.0185616) q[2];
sx q[2];
rz(-1.8862714) q[2];
sx q[2];
rz(-1.4184784) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.84049837) q[1];
sx q[1];
rz(-2.4972895) q[1];
sx q[1];
rz(-1.3703764) q[1];
x q[2];
rz(-2.850344) q[3];
sx q[3];
rz(-2.2403803) q[3];
sx q[3];
rz(2.254829) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.5608998) q[2];
sx q[2];
rz(-2.7036724) q[2];
sx q[2];
rz(1.8221347) q[2];
rz(-0.27160078) q[3];
sx q[3];
rz(-1.8235794) q[3];
sx q[3];
rz(1.1118838) q[3];
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
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.44529799) q[0];
sx q[0];
rz(-2.2309208) q[0];
sx q[0];
rz(-1.4060422) q[0];
rz(1.4259074) q[1];
sx q[1];
rz(-2.0366663) q[1];
sx q[1];
rz(-1.8941194) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46215993) q[0];
sx q[0];
rz(-1.3094256) q[0];
sx q[0];
rz(-2.0716459) q[0];
rz(-pi) q[1];
x q[1];
rz(1.6883759) q[2];
sx q[2];
rz(-1.5167859) q[2];
sx q[2];
rz(-2.5117995) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.61037072) q[1];
sx q[1];
rz(-2.7970813) q[1];
sx q[1];
rz(0.12126874) q[1];
rz(-pi) q[2];
rz(-2.3214809) q[3];
sx q[3];
rz(-0.81104987) q[3];
sx q[3];
rz(0.91203989) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.4282816) q[2];
sx q[2];
rz(-0.2321299) q[2];
sx q[2];
rz(3.0408119) q[2];
rz(3.1320324) q[3];
sx q[3];
rz(-1.5567501) q[3];
sx q[3];
rz(-1.6077707) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(2.1495001) q[0];
sx q[0];
rz(-1.7763573) q[0];
sx q[0];
rz(-1.1471164) q[0];
rz(-2.5485582) q[1];
sx q[1];
rz(-1.7094163) q[1];
sx q[1];
rz(-0.97547668) q[1];
rz(1.7375853) q[2];
sx q[2];
rz(-0.40344147) q[2];
sx q[2];
rz(-0.16735195) q[2];
rz(1.5884052) q[3];
sx q[3];
rz(-1.3497769) q[3];
sx q[3];
rz(-2.6365437) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
