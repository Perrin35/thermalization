OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.0176528) q[0];
sx q[0];
rz(-2.2247563) q[0];
sx q[0];
rz(0.43489781) q[0];
rz(-1.544156) q[1];
sx q[1];
rz(-0.51900744) q[1];
sx q[1];
rz(-2.5595698) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3618185) q[0];
sx q[0];
rz(-1.7248745) q[0];
sx q[0];
rz(1.3067354) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.97442128) q[2];
sx q[2];
rz(-2.4108464) q[2];
sx q[2];
rz(1.3517018) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.37108818) q[1];
sx q[1];
rz(-0.50482115) q[1];
sx q[1];
rz(-2.0981028) q[1];
rz(-pi) q[2];
rz(-0.1780171) q[3];
sx q[3];
rz(-1.8478763) q[3];
sx q[3];
rz(-1.2704364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.5883098) q[2];
sx q[2];
rz(-1.2574235) q[2];
sx q[2];
rz(2.8498939) q[2];
rz(-0.45082539) q[3];
sx q[3];
rz(-2.8901633) q[3];
sx q[3];
rz(-2.5341471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70158231) q[0];
sx q[0];
rz(-0.99047438) q[0];
sx q[0];
rz(1.095358) q[0];
rz(1.9502684) q[1];
sx q[1];
rz(-2.2215863) q[1];
sx q[1];
rz(2.2390168) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13690925) q[0];
sx q[0];
rz(-2.1666514) q[0];
sx q[0];
rz(2.3803902) q[0];
rz(-pi) q[1];
rz(-1.9368725) q[2];
sx q[2];
rz(-0.70209018) q[2];
sx q[2];
rz(3.0099208) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.7831259) q[1];
sx q[1];
rz(-1.2791113) q[1];
sx q[1];
rz(-1.6963708) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.91530494) q[3];
sx q[3];
rz(-2.4065131) q[3];
sx q[3];
rz(-0.32441586) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-0.1878745) q[2];
sx q[2];
rz(-0.97430054) q[2];
sx q[2];
rz(0.11889674) q[2];
rz(-2.4925354) q[3];
sx q[3];
rz(-0.18638149) q[3];
sx q[3];
rz(1.4547179) q[3];
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
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4897937) q[0];
sx q[0];
rz(-1.8875903) q[0];
sx q[0];
rz(0.5157665) q[0];
rz(-2.0491397) q[1];
sx q[1];
rz(-2.7383883) q[1];
sx q[1];
rz(-1.2752424) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9495205) q[0];
sx q[0];
rz(-1.6123591) q[0];
sx q[0];
rz(1.5526718) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.8280806) q[2];
sx q[2];
rz(-1.7898774) q[2];
sx q[2];
rz(0.30922019) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0787492) q[1];
sx q[1];
rz(-2.242452) q[1];
sx q[1];
rz(0.25154227) q[1];
rz(-pi) q[2];
rz(-3.0413275) q[3];
sx q[3];
rz(-0.79660049) q[3];
sx q[3];
rz(-1.0880566) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.0755997) q[2];
sx q[2];
rz(-1.3174572) q[2];
sx q[2];
rz(-1.1598178) q[2];
rz(-2.6521111) q[3];
sx q[3];
rz(-1.5882086) q[3];
sx q[3];
rz(2.421853) q[3];
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
rz(-pi/2) q[0];
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2751665) q[0];
sx q[0];
rz(-2.8186099) q[0];
sx q[0];
rz(0.48165709) q[0];
rz(-2.2109168) q[1];
sx q[1];
rz(-1.8447256) q[1];
sx q[1];
rz(0.01893386) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.16068072) q[0];
sx q[0];
rz(-1.5917516) q[0];
sx q[0];
rz(1.5207855) q[0];
rz(-pi) q[1];
rz(-0.04425116) q[2];
sx q[2];
rz(-2.2328064) q[2];
sx q[2];
rz(0.64427441) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.1065797) q[1];
sx q[1];
rz(-1.0570608) q[1];
sx q[1];
rz(-0.68961838) q[1];
rz(-pi) q[2];
rz(2.3036257) q[3];
sx q[3];
rz(-0.95221704) q[3];
sx q[3];
rz(-3.0493506) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.9229752) q[2];
sx q[2];
rz(-0.3564035) q[2];
sx q[2];
rz(-0.74529988) q[2];
rz(-0.37799147) q[3];
sx q[3];
rz(-1.9972921) q[3];
sx q[3];
rz(2.2530344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.22076631) q[0];
sx q[0];
rz(-2.8296318) q[0];
sx q[0];
rz(-1.1908603) q[0];
rz(2.6530755) q[1];
sx q[1];
rz(-0.91594511) q[1];
sx q[1];
rz(0.8078422) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61208439) q[0];
sx q[0];
rz(-2.8017462) q[0];
sx q[0];
rz(-1.8456949) q[0];
x q[1];
rz(2.0338397) q[2];
sx q[2];
rz(-1.7757738) q[2];
sx q[2];
rz(-3.0489655) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.1792308) q[1];
sx q[1];
rz(-0.46474248) q[1];
sx q[1];
rz(2.9443113) q[1];
rz(-pi) q[2];
x q[2];
rz(2.1937815) q[3];
sx q[3];
rz(-1.4075507) q[3];
sx q[3];
rz(2.3415274) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.1390344) q[2];
sx q[2];
rz(-0.99271861) q[2];
sx q[2];
rz(0.16415088) q[2];
rz(-0.74603355) q[3];
sx q[3];
rz(-1.7697216) q[3];
sx q[3];
rz(0.073808864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98110759) q[0];
sx q[0];
rz(-1.8758513) q[0];
sx q[0];
rz(-0.72738457) q[0];
rz(-0.47850594) q[1];
sx q[1];
rz(-1.6886657) q[1];
sx q[1];
rz(-0.7368288) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3664624) q[0];
sx q[0];
rz(-1.7084165) q[0];
sx q[0];
rz(-0.055995106) q[0];
rz(-pi) q[1];
rz(1.9987891) q[2];
sx q[2];
rz(-2.8341132) q[2];
sx q[2];
rz(-0.82915074) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.9068272) q[1];
sx q[1];
rz(-1.9363469) q[1];
sx q[1];
rz(0.91616456) q[1];
rz(-pi) q[2];
x q[2];
rz(3.1276631) q[3];
sx q[3];
rz(-1.1831814) q[3];
sx q[3];
rz(0.25957169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.14787802) q[2];
sx q[2];
rz(-1.0241877) q[2];
sx q[2];
rz(-0.68515879) q[2];
rz(2.1327175) q[3];
sx q[3];
rz(-0.69005552) q[3];
sx q[3];
rz(-2.8907997) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.93609) q[0];
sx q[0];
rz(-3.0476373) q[0];
sx q[0];
rz(2.0482546) q[0];
rz(3.1178442) q[1];
sx q[1];
rz(-0.7370342) q[1];
sx q[1];
rz(2.645983) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0422875) q[0];
sx q[0];
rz(-1.6368425) q[0];
sx q[0];
rz(1.582986) q[0];
x q[1];
rz(-2.7205207) q[2];
sx q[2];
rz(-0.47343987) q[2];
sx q[2];
rz(0.978549) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-1.1417947) q[1];
sx q[1];
rz(-0.3965946) q[1];
sx q[1];
rz(-0.48626089) q[1];
x q[2];
rz(2.2836693) q[3];
sx q[3];
rz(-0.61238785) q[3];
sx q[3];
rz(-0.96549406) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2451943) q[2];
sx q[2];
rz(-2.3485025) q[2];
sx q[2];
rz(0.29655656) q[2];
rz(-0.5101997) q[3];
sx q[3];
rz(-1.5254131) q[3];
sx q[3];
rz(2.8701674) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9940014) q[0];
sx q[0];
rz(-2.1864102) q[0];
sx q[0];
rz(-1.1543132) q[0];
rz(1.9860024) q[1];
sx q[1];
rz(-2.4842333) q[1];
sx q[1];
rz(-1.4591699) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4765035) q[0];
sx q[0];
rz(-0.71528331) q[0];
sx q[0];
rz(-3.1066549) q[0];
rz(-pi) q[1];
rz(1.5695249) q[2];
sx q[2];
rz(-2.4172815) q[2];
sx q[2];
rz(1.5226165) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.1016804) q[1];
sx q[1];
rz(-1.7500568) q[1];
sx q[1];
rz(2.0800566) q[1];
rz(-pi) q[2];
rz(-2.5063839) q[3];
sx q[3];
rz(-2.016474) q[3];
sx q[3];
rz(-3.095568) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3079754) q[2];
sx q[2];
rz(-0.54215446) q[2];
sx q[2];
rz(1.3758434) q[2];
rz(-3.0552676) q[3];
sx q[3];
rz(-2.3089246) q[3];
sx q[3];
rz(1.255792) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
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
rz(-1.9432705) q[0];
sx q[0];
rz(-2.4477796) q[0];
sx q[0];
rz(-0.86945239) q[0];
rz(2.4122639) q[1];
sx q[1];
rz(-0.84356934) q[1];
sx q[1];
rz(0.23390153) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9651523) q[0];
sx q[0];
rz(-1.5230706) q[0];
sx q[0];
rz(-0.41712572) q[0];
rz(-0.42371427) q[2];
sx q[2];
rz(-0.22988453) q[2];
sx q[2];
rz(-0.73424852) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.42185509) q[1];
sx q[1];
rz(-0.46793391) q[1];
sx q[1];
rz(0.33166064) q[1];
rz(-1.2083268) q[3];
sx q[3];
rz(-0.7504645) q[3];
sx q[3];
rz(-3.124231) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(3.0539315) q[2];
sx q[2];
rz(-1.4940741) q[2];
sx q[2];
rz(1.9653448) q[2];
rz(0.67115274) q[3];
sx q[3];
rz(-1.2495557) q[3];
sx q[3];
rz(1.6254856) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.43435708) q[0];
sx q[0];
rz(-0.86152995) q[0];
sx q[0];
rz(1.0435411) q[0];
rz(-0.097298233) q[1];
sx q[1];
rz(-2.0847335) q[1];
sx q[1];
rz(1.1204488) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1277043) q[0];
sx q[0];
rz(-0.59895016) q[0];
sx q[0];
rz(1.6779792) q[0];
rz(-pi) q[1];
rz(-2.3515969) q[2];
sx q[2];
rz(-0.95271356) q[2];
sx q[2];
rz(0.8481942) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.7880558) q[1];
sx q[1];
rz(-2.9573943) q[1];
sx q[1];
rz(2.8723529) q[1];
rz(-pi) q[2];
rz(2.9531526) q[3];
sx q[3];
rz(-2.0456637) q[3];
sx q[3];
rz(1.4120073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.4572767) q[2];
sx q[2];
rz(-2.5280759) q[2];
sx q[2];
rz(1.3555917) q[2];
rz(-1.5654303) q[3];
sx q[3];
rz(-1.0836982) q[3];
sx q[3];
rz(2.1255597) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.23715699) q[0];
sx q[0];
rz(-0.56023993) q[0];
sx q[0];
rz(2.3518363) q[0];
rz(0.65658983) q[1];
sx q[1];
rz(-1.1142535) q[1];
sx q[1];
rz(2.5744892) q[1];
rz(-0.97548998) q[2];
sx q[2];
rz(-2.3885623) q[2];
sx q[2];
rz(1.9293712) q[2];
rz(2.5593467) q[3];
sx q[3];
rz(-2.7979294) q[3];
sx q[3];
rz(2.9002849) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
