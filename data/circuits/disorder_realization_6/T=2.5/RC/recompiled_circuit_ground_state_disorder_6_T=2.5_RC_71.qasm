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
rz(4.1832357) q[1];
sx q[1];
rz(15.19419) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0002961) q[0];
sx q[0];
rz(-1.0736672) q[0];
sx q[0];
rz(0.5750467) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.47494048) q[2];
sx q[2];
rz(-1.5154334) q[2];
sx q[2];
rz(-0.1794912) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.7440255) q[1];
sx q[1];
rz(-1.7850707) q[1];
sx q[1];
rz(2.812723) q[1];
rz(-pi) q[2];
rz(-0.051187201) q[3];
sx q[3];
rz(-0.59268206) q[3];
sx q[3];
rz(0.41646233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.4673246) q[2];
sx q[2];
rz(-1.5200204) q[2];
sx q[2];
rz(-2.8094214) q[2];
rz(2.8915571) q[3];
sx q[3];
rz(-1.9146405) q[3];
sx q[3];
rz(1.2927607) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.2501204) q[0];
sx q[0];
rz(-1.253506) q[0];
sx q[0];
rz(-0.44723311) q[0];
rz(-2.1121292) q[1];
sx q[1];
rz(-2.5599458) q[1];
sx q[1];
rz(-0.10118016) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6529236) q[0];
sx q[0];
rz(-2.5258668) q[0];
sx q[0];
rz(0.73955817) q[0];
x q[1];
rz(1.871505) q[2];
sx q[2];
rz(-2.5413587) q[2];
sx q[2];
rz(2.9832961) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.78522462) q[1];
sx q[1];
rz(-1.1572946) q[1];
sx q[1];
rz(1.2797536) q[1];
rz(-pi) q[2];
rz(1.7216136) q[3];
sx q[3];
rz(-2.4972417) q[3];
sx q[3];
rz(2.2805813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.0589361) q[2];
sx q[2];
rz(-0.75581789) q[2];
sx q[2];
rz(1.6015046) q[2];
rz(-0.61795175) q[3];
sx q[3];
rz(-2.5585585) q[3];
sx q[3];
rz(1.9790953) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.59442941) q[0];
sx q[0];
rz(-1.0048486) q[0];
sx q[0];
rz(-0.51602236) q[0];
rz(-2.0891321) q[1];
sx q[1];
rz(-2.2155589) q[1];
sx q[1];
rz(1.9693536) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2207551) q[0];
sx q[0];
rz(-1.1467548) q[0];
sx q[0];
rz(1.5586669) q[0];
rz(3.0686106) q[2];
sx q[2];
rz(-0.86162815) q[2];
sx q[2];
rz(0.87984171) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7683923) q[1];
sx q[1];
rz(-0.52881634) q[1];
sx q[1];
rz(1.1952728) q[1];
x q[2];
rz(0.013962176) q[3];
sx q[3];
rz(-0.80028557) q[3];
sx q[3];
rz(1.6228907) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8664794) q[2];
sx q[2];
rz(-1.8786414) q[2];
sx q[2];
rz(2.4647554) q[2];
rz(0.46946851) q[3];
sx q[3];
rz(-0.89478409) q[3];
sx q[3];
rz(1.9278056) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6770099) q[0];
sx q[0];
rz(-1.0053758) q[0];
sx q[0];
rz(1.9418035) q[0];
rz(2.5489573) q[1];
sx q[1];
rz(-2.0694144) q[1];
sx q[1];
rz(1.6823654) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1593247) q[0];
sx q[0];
rz(-1.3749116) q[0];
sx q[0];
rz(2.5513493) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8921669) q[2];
sx q[2];
rz(-2.327033) q[2];
sx q[2];
rz(0.047175353) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.063733405) q[1];
sx q[1];
rz(-0.95653906) q[1];
sx q[1];
rz(0.31645223) q[1];
x q[2];
rz(-2.4577971) q[3];
sx q[3];
rz(-1.8701564) q[3];
sx q[3];
rz(-1.2049068) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
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
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
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
rz(0.84363627) q[0];
sx q[0];
rz(-1.4763259) q[0];
sx q[0];
rz(0.27808878) q[0];
rz(0.9043215) q[1];
sx q[1];
rz(-1.6878637) q[1];
sx q[1];
rz(-2.1542737) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.41288677) q[0];
sx q[0];
rz(-2.1653919) q[0];
sx q[0];
rz(-2.8534858) q[0];
rz(1.5989207) q[2];
sx q[2];
rz(-1.7076945) q[2];
sx q[2];
rz(1.9134476) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6254355) q[1];
sx q[1];
rz(-1.9599116) q[1];
sx q[1];
rz(0.066246943) q[1];
rz(-pi) q[2];
x q[2];
rz(1.1313391) q[3];
sx q[3];
rz(-2.5601697) q[3];
sx q[3];
rz(-1.9245337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.014293369) q[2];
sx q[2];
rz(-0.96724808) q[2];
sx q[2];
rz(-2.3822752) q[2];
rz(1.6589818) q[3];
sx q[3];
rz(-1.9448152) q[3];
sx q[3];
rz(2.7827941) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9229014) q[0];
sx q[0];
rz(-2.3221115) q[0];
sx q[0];
rz(-0.6915834) q[0];
rz(-2.2870731) q[1];
sx q[1];
rz(-0.71047345) q[1];
sx q[1];
rz(2.7202594) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79165073) q[0];
sx q[0];
rz(-1.5022819) q[0];
sx q[0];
rz(-1.1380026) q[0];
rz(1.2796655) q[2];
sx q[2];
rz(-1.8270681) q[2];
sx q[2];
rz(2.1879183) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.2273352) q[1];
sx q[1];
rz(-0.37082878) q[1];
sx q[1];
rz(-0.23359681) q[1];
rz(-pi) q[2];
rz(-2.5365127) q[3];
sx q[3];
rz(-1.7972094) q[3];
sx q[3];
rz(0.070158557) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.93115807) q[2];
sx q[2];
rz(-1.4782108) q[2];
sx q[2];
rz(0.20273905) q[2];
rz(1.5431131) q[3];
sx q[3];
rz(-0.32035443) q[3];
sx q[3];
rz(1.4690442) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8572674) q[0];
sx q[0];
rz(-1.7131282) q[0];
sx q[0];
rz(2.7427234) q[0];
rz(1.9786037) q[1];
sx q[1];
rz(-2.3154924) q[1];
sx q[1];
rz(-0.2805447) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3464081) q[0];
sx q[0];
rz(-1.9530256) q[0];
sx q[0];
rz(0.43117492) q[0];
rz(-pi) q[1];
x q[1];
rz(0.89338867) q[2];
sx q[2];
rz(-1.9147493) q[2];
sx q[2];
rz(2.266573) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.0979586) q[1];
sx q[1];
rz(-0.53501883) q[1];
sx q[1];
rz(-0.94162264) q[1];
rz(-pi) q[2];
rz(-0.97784247) q[3];
sx q[3];
rz(-2.2740318) q[3];
sx q[3];
rz(-0.11444005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.82105381) q[2];
sx q[2];
rz(-0.97964779) q[2];
sx q[2];
rz(-1.6542356) q[2];
rz(-2.1982543) q[3];
sx q[3];
rz(-2.2127559) q[3];
sx q[3];
rz(-1.1613891) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5281552) q[0];
sx q[0];
rz(-1.4035839) q[0];
sx q[0];
rz(-0.50233895) q[0];
rz(-2.595064) q[1];
sx q[1];
rz(-1.2556475) q[1];
sx q[1];
rz(-0.46195236) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.079507) q[0];
sx q[0];
rz(-1.8621729) q[0];
sx q[0];
rz(2.0847005) q[0];
x q[1];
rz(-1.6335604) q[2];
sx q[2];
rz(-2.4893005) q[2];
sx q[2];
rz(-1.4029274) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.6071872) q[1];
sx q[1];
rz(-2.7173373) q[1];
sx q[1];
rz(-1.8380182) q[1];
x q[2];
rz(1.5580719) q[3];
sx q[3];
rz(-2.1012348) q[3];
sx q[3];
rz(-2.1761571) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.168557) q[2];
sx q[2];
rz(-2.904197) q[2];
sx q[2];
rz(2.9711704) q[2];
rz(2.6863344) q[3];
sx q[3];
rz(-1.596343) q[3];
sx q[3];
rz(-2.4236603) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
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
rz(0.43389854) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.949373) q[0];
sx q[0];
rz(-1.1447971) q[0];
sx q[0];
rz(0.98652391) q[0];
x q[1];
rz(-0.92484464) q[2];
sx q[2];
rz(-2.6000823) q[2];
sx q[2];
rz(-0.72586593) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.89131195) q[1];
sx q[1];
rz(-1.4509333) q[1];
sx q[1];
rz(-0.93617546) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.9187886) q[3];
sx q[3];
rz(-0.7211313) q[3];
sx q[3];
rz(-2.7047771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.58069289) q[2];
sx q[2];
rz(-0.43792024) q[2];
sx q[2];
rz(-1.8221347) q[2];
rz(-2.8699919) q[3];
sx q[3];
rz(-1.3180132) q[3];
sx q[3];
rz(1.1118838) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi/2) q[3];
sx q[3];
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
rz(-2.6962947) q[0];
sx q[0];
rz(-2.2309208) q[0];
sx q[0];
rz(-1.7355504) q[0];
rz(1.7156853) q[1];
sx q[1];
rz(-1.1049263) q[1];
sx q[1];
rz(-1.8941194) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5496779) q[0];
sx q[0];
rz(-2.5818338) q[0];
sx q[0];
rz(1.0625398) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.054385238) q[2];
sx q[2];
rz(-1.6882036) q[2];
sx q[2];
rz(0.94737999) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0669603) q[1];
sx q[1];
rz(-1.5299284) q[1];
sx q[1];
rz(0.34217477) q[1];
rz(-0.62274751) q[3];
sx q[3];
rz(-2.1295665) q[3];
sx q[3];
rz(-1.2947242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.4282816) q[2];
sx q[2];
rz(-0.2321299) q[2];
sx q[2];
rz(-3.0408119) q[2];
rz(-0.0095602592) q[3];
sx q[3];
rz(-1.5848426) q[3];
sx q[3];
rz(1.6077707) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.99209256) q[0];
sx q[0];
rz(-1.7763573) q[0];
sx q[0];
rz(-1.1471164) q[0];
rz(-0.59303444) q[1];
sx q[1];
rz(-1.4321764) q[1];
sx q[1];
rz(2.166116) q[1];
rz(-1.4040074) q[2];
sx q[2];
rz(-0.40344147) q[2];
sx q[2];
rz(-0.16735195) q[2];
rz(0.22105263) q[3];
sx q[3];
rz(-1.5879769) q[3];
sx q[3];
rz(-1.069608) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
