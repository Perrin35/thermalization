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
rz(6.8721539) q[0];
sx q[0];
rz(6.4018259) q[0];
rz(2.150382) q[1];
sx q[1];
rz(-1.041643) q[1];
sx q[1];
rz(-2.6278194) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27072752) q[0];
sx q[0];
rz(-1.0723734) q[0];
sx q[0];
rz(0.99683572) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0209885) q[2];
sx q[2];
rz(-0.47791156) q[2];
sx q[2];
rz(-1.284104) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.5255907) q[1];
sx q[1];
rz(-2.7512128) q[1];
sx q[1];
rz(2.5486773) q[1];
rz(-2.5495177) q[3];
sx q[3];
rz(-1.5993803) q[3];
sx q[3];
rz(2.0297272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.674268) q[2];
sx q[2];
rz(-1.5200204) q[2];
sx q[2];
rz(2.8094214) q[2];
rz(2.8915571) q[3];
sx q[3];
rz(-1.2269521) q[3];
sx q[3];
rz(1.8488319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8914723) q[0];
sx q[0];
rz(-1.8880867) q[0];
sx q[0];
rz(-2.6943595) q[0];
rz(2.1121292) q[1];
sx q[1];
rz(-0.58164683) q[1];
sx q[1];
rz(-0.10118016) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6995938) q[0];
sx q[0];
rz(-1.170982) q[0];
sx q[0];
rz(-2.6599822) q[0];
x q[1];
rz(2.9415628) q[2];
sx q[2];
rz(-1.0009655) q[2];
sx q[2];
rz(-0.20117682) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.78522462) q[1];
sx q[1];
rz(-1.984298) q[1];
sx q[1];
rz(-1.861839) q[1];
rz(-pi) q[2];
rz(-0.93191913) q[3];
sx q[3];
rz(-1.4804236) q[3];
sx q[3];
rz(-0.83069176) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.0589361) q[2];
sx q[2];
rz(-0.75581789) q[2];
sx q[2];
rz(1.6015046) q[2];
rz(-2.5236409) q[3];
sx q[3];
rz(-0.58303419) q[3];
sx q[3];
rz(1.9790953) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(0.59442941) q[0];
sx q[0];
rz(-2.1367441) q[0];
sx q[0];
rz(0.51602236) q[0];
rz(-1.0524606) q[1];
sx q[1];
rz(-0.92603374) q[1];
sx q[1];
rz(1.9693536) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7965423) q[0];
sx q[0];
rz(-1.5818514) q[0];
sx q[0];
rz(-2.7175236) q[0];
x q[1];
rz(0.86031057) q[2];
sx q[2];
rz(-1.5154308) q[2];
sx q[2];
rz(-2.4982128) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7683923) q[1];
sx q[1];
rz(-2.6127763) q[1];
sx q[1];
rz(-1.1952728) q[1];
x q[2];
rz(1.5564136) q[3];
sx q[3];
rz(-2.3709815) q[3];
sx q[3];
rz(1.5387467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8664794) q[2];
sx q[2];
rz(-1.8786414) q[2];
sx q[2];
rz(-0.67683721) q[2];
rz(2.6721241) q[3];
sx q[3];
rz(-0.89478409) q[3];
sx q[3];
rz(-1.9278056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4645828) q[0];
sx q[0];
rz(-1.0053758) q[0];
sx q[0];
rz(-1.1997892) q[0];
rz(0.59263539) q[1];
sx q[1];
rz(-2.0694144) q[1];
sx q[1];
rz(1.4592272) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1593247) q[0];
sx q[0];
rz(-1.3749116) q[0];
sx q[0];
rz(0.59024337) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3427354) q[2];
sx q[2];
rz(-1.3902544) q[2];
sx q[2];
rz(-1.7910166) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.063733405) q[1];
sx q[1];
rz(-2.1850536) q[1];
sx q[1];
rz(2.8251404) q[1];
rz(-pi) q[2];
rz(-0.68379559) q[3];
sx q[3];
rz(-1.8701564) q[3];
sx q[3];
rz(-1.9366858) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.6478641) q[2];
sx q[2];
rz(-2.04546) q[2];
sx q[2];
rz(0.68191051) q[2];
rz(2.72825) q[3];
sx q[3];
rz(-1.5324493) q[3];
sx q[3];
rz(-1.9857064) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2979564) q[0];
sx q[0];
rz(-1.6652668) q[0];
sx q[0];
rz(0.27808878) q[0];
rz(0.9043215) q[1];
sx q[1];
rz(-1.6878637) q[1];
sx q[1];
rz(-2.1542737) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7287059) q[0];
sx q[0];
rz(-2.1653919) q[0];
sx q[0];
rz(0.28810687) q[0];
rz(-0.1369517) q[2];
sx q[2];
rz(-1.5429351) q[2];
sx q[2];
rz(2.7951023) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(1.6254355) q[1];
sx q[1];
rz(-1.9599116) q[1];
sx q[1];
rz(-0.066246943) q[1];
rz(-pi) q[2];
rz(1.0342399) q[3];
sx q[3];
rz(-1.3349541) q[3];
sx q[3];
rz(0.72808108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.014293369) q[2];
sx q[2];
rz(-0.96724808) q[2];
sx q[2];
rz(2.3822752) q[2];
rz(-1.6589818) q[3];
sx q[3];
rz(-1.1967775) q[3];
sx q[3];
rz(-0.35879859) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.21869126) q[0];
sx q[0];
rz(-0.8194812) q[0];
sx q[0];
rz(-0.6915834) q[0];
rz(2.2870731) q[1];
sx q[1];
rz(-0.71047345) q[1];
sx q[1];
rz(0.42133322) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2153361) q[0];
sx q[0];
rz(-2.7037482) q[0];
sx q[0];
rz(-1.4086176) q[0];
rz(-pi) q[1];
rz(2.3106903) q[2];
sx q[2];
rz(-0.38543265) q[2];
sx q[2];
rz(1.3192434) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(3.1189733) q[1];
sx q[1];
rz(-1.9310863) q[1];
sx q[1];
rz(-1.4810356) q[1];
rz(-2.756802) q[3];
sx q[3];
rz(-0.64105415) q[3];
sx q[3];
rz(1.8144864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2104346) q[2];
sx q[2];
rz(-1.6633818) q[2];
sx q[2];
rz(0.20273905) q[2];
rz(1.5984795) q[3];
sx q[3];
rz(-2.8212382) q[3];
sx q[3];
rz(1.4690442) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8572674) q[0];
sx q[0];
rz(-1.4284644) q[0];
sx q[0];
rz(-2.7427234) q[0];
rz(-1.162989) q[1];
sx q[1];
rz(-2.3154924) q[1];
sx q[1];
rz(-0.2805447) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4568854) q[0];
sx q[0];
rz(-2.5734757) q[0];
sx q[0];
rz(2.3756308) q[0];
rz(-pi) q[1];
rz(2.7107096) q[2];
sx q[2];
rz(-0.93965778) q[2];
sx q[2];
rz(0.43090303) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.1093367) q[1];
sx q[1];
rz(-1.2660626) q[1];
sx q[1];
rz(-2.017657) q[1];
rz(-pi) q[2];
rz(0.97784247) q[3];
sx q[3];
rz(-2.2740318) q[3];
sx q[3];
rz(0.11444005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.82105381) q[2];
sx q[2];
rz(-0.97964779) q[2];
sx q[2];
rz(1.6542356) q[2];
rz(-0.94333831) q[3];
sx q[3];
rz(-2.2127559) q[3];
sx q[3];
rz(-1.9802035) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.6134375) q[0];
sx q[0];
rz(-1.4035839) q[0];
sx q[0];
rz(-0.50233895) q[0];
rz(-2.595064) q[1];
sx q[1];
rz(-1.2556475) q[1];
sx q[1];
rz(2.6796403) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7936458) q[0];
sx q[0];
rz(-1.0805305) q[0];
sx q[0];
rz(-2.8099174) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0937198) q[2];
sx q[2];
rz(-2.2215853) q[2];
sx q[2];
rz(1.4818459) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.71910243) q[1];
sx q[1];
rz(-1.6797069) q[1];
sx q[1];
rz(-1.1599354) q[1];
rz(-pi) q[2];
x q[2];
rz(0.53047385) q[3];
sx q[3];
rz(-1.5817721) q[3];
sx q[3];
rz(2.5297942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.168557) q[2];
sx q[2];
rz(-2.904197) q[2];
sx q[2];
rz(-2.9711704) q[2];
rz(-0.45525822) q[3];
sx q[3];
rz(-1.596343) q[3];
sx q[3];
rz(0.71793238) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.97487226) q[0];
sx q[0];
rz(-1.6396739) q[0];
sx q[0];
rz(-2.9611452) q[0];
rz(0.51668733) q[1];
sx q[1];
rz(-1.5770715) q[1];
sx q[1];
rz(-2.7076941) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2045) q[0];
sx q[0];
rz(-0.70815101) q[0];
sx q[0];
rz(-2.2592179) q[0];
rz(-pi) q[1];
rz(-2.0185616) q[2];
sx q[2];
rz(-1.2553213) q[2];
sx q[2];
rz(-1.7231143) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.3010943) q[1];
sx q[1];
rz(-0.64430311) q[1];
sx q[1];
rz(1.3703764) q[1];
rz(2.2614165) q[3];
sx q[3];
rz(-1.7978906) q[3];
sx q[3];
rz(-0.86798641) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5608998) q[2];
sx q[2];
rz(-2.7036724) q[2];
sx q[2];
rz(1.3194579) q[2];
rz(-2.8699919) q[3];
sx q[3];
rz(-1.3180132) q[3];
sx q[3];
rz(-2.0297089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.44529799) q[0];
sx q[0];
rz(-2.2309208) q[0];
sx q[0];
rz(-1.7355504) q[0];
rz(-1.7156853) q[1];
sx q[1];
rz(-2.0366663) q[1];
sx q[1];
rz(-1.8941194) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5919148) q[0];
sx q[0];
rz(-0.55975883) q[0];
sx q[0];
rz(-1.0625398) q[0];
rz(-pi) q[1];
rz(-2.0026454) q[2];
sx q[2];
rz(-3.0122535) q[2];
sx q[2];
rz(-0.51233385) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.48161246) q[1];
sx q[1];
rz(-1.9126737) q[1];
sx q[1];
rz(1.6141763) q[1];
x q[2];
rz(0.62274751) q[3];
sx q[3];
rz(-2.1295665) q[3];
sx q[3];
rz(-1.8468685) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.71331104) q[2];
sx q[2];
rz(-2.9094628) q[2];
sx q[2];
rz(0.10078079) q[2];
rz(-3.1320324) q[3];
sx q[3];
rz(-1.5567501) q[3];
sx q[3];
rz(-1.5338219) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1495001) q[0];
sx q[0];
rz(-1.3652353) q[0];
sx q[0];
rz(1.9944763) q[0];
rz(-2.5485582) q[1];
sx q[1];
rz(-1.7094163) q[1];
sx q[1];
rz(-0.97547668) q[1];
rz(-3.0708457) q[2];
sx q[2];
rz(-1.9683217) q[2];
sx q[2];
rz(-0.34839658) q[2];
rz(-1.5884052) q[3];
sx q[3];
rz(-1.7918158) q[3];
sx q[3];
rz(0.50504897) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
