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
rz(2.0496378) q[0];
sx q[0];
rz(-2.0593934) q[0];
sx q[0];
rz(1.9424633) q[0];
rz(-2.8973329) q[1];
sx q[1];
rz(-1.8241939) q[1];
sx q[1];
rz(0.86950818) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8438709) q[0];
sx q[0];
rz(-1.3605357) q[0];
sx q[0];
rz(-0.79550996) q[0];
rz(-pi) q[1];
rz(2.7156285) q[2];
sx q[2];
rz(-1.1168343) q[2];
sx q[2];
rz(-0.23655836) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(3.0119734) q[1];
sx q[1];
rz(-1.8070613) q[1];
sx q[1];
rz(-0.70567997) q[1];
rz(0.77750364) q[3];
sx q[3];
rz(-2.1453259) q[3];
sx q[3];
rz(-0.053701775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.96183744) q[2];
sx q[2];
rz(-1.1575674) q[2];
sx q[2];
rz(-1.301514) q[2];
rz(-2.2573722) q[3];
sx q[3];
rz(-1.0853465) q[3];
sx q[3];
rz(1.4475383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
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
rz(-2.7262511) q[0];
sx q[0];
rz(-1.2301507) q[0];
sx q[0];
rz(-3.0905261) q[0];
rz(0.49634936) q[1];
sx q[1];
rz(-1.2800848) q[1];
sx q[1];
rz(3.0275717) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89897777) q[0];
sx q[0];
rz(-2.0847843) q[0];
sx q[0];
rz(1.32715) q[0];
rz(-pi) q[1];
x q[1];
rz(0.013234303) q[2];
sx q[2];
rz(-0.87761384) q[2];
sx q[2];
rz(-0.40738328) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.025578827) q[1];
sx q[1];
rz(-2.4878628) q[1];
sx q[1];
rz(0.81029256) q[1];
x q[2];
rz(2.5258371) q[3];
sx q[3];
rz(-2.1020707) q[3];
sx q[3];
rz(2.5221276) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.8337635) q[2];
sx q[2];
rz(-1.2276063) q[2];
sx q[2];
rz(0.366079) q[2];
rz(2.6858373) q[3];
sx q[3];
rz(-2.3594806) q[3];
sx q[3];
rz(2.8428049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.61641055) q[0];
sx q[0];
rz(-0.95935482) q[0];
sx q[0];
rz(2.8042703) q[0];
rz(1.6711383) q[1];
sx q[1];
rz(-2.2098139) q[1];
sx q[1];
rz(-3.0477188) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.84067837) q[0];
sx q[0];
rz(-0.72193679) q[0];
sx q[0];
rz(3.0283699) q[0];
rz(-pi) q[1];
rz(-2.6930935) q[2];
sx q[2];
rz(-1.8828159) q[2];
sx q[2];
rz(-1.2879368) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.17901668) q[1];
sx q[1];
rz(-1.1454795) q[1];
sx q[1];
rz(-3.1326181) q[1];
rz(-2.8416614) q[3];
sx q[3];
rz(-0.75977221) q[3];
sx q[3];
rz(-2.531736) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.9570534) q[2];
sx q[2];
rz(-0.55593714) q[2];
sx q[2];
rz(-2.8957193) q[2];
rz(-2.9546514) q[3];
sx q[3];
rz(-1.4176466) q[3];
sx q[3];
rz(2.7174182) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
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
rz(-2.4251855) q[0];
sx q[0];
rz(-2.462429) q[0];
sx q[0];
rz(2.9492522) q[0];
rz(1.1771857) q[1];
sx q[1];
rz(-2.3277551) q[1];
sx q[1];
rz(0.40843931) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.724736) q[0];
sx q[0];
rz(-1.4166066) q[0];
sx q[0];
rz(1.6282999) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.8544418) q[2];
sx q[2];
rz(-0.39098323) q[2];
sx q[2];
rz(1.4650844) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7725984) q[1];
sx q[1];
rz(-2.1647758) q[1];
sx q[1];
rz(1.0812378) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.2564234) q[3];
sx q[3];
rz(-1.4082296) q[3];
sx q[3];
rz(-3.0522961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.283215) q[2];
sx q[2];
rz(-2.242531) q[2];
sx q[2];
rz(-2.9646207) q[2];
rz(-2.5519798) q[3];
sx q[3];
rz(-1.4859345) q[3];
sx q[3];
rz(-1.0346712) q[3];
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
sx q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4347587) q[0];
sx q[0];
rz(-0.85586268) q[0];
sx q[0];
rz(2.6245497) q[0];
rz(0.21273461) q[1];
sx q[1];
rz(-2.0566302) q[1];
sx q[1];
rz(0.83232602) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34236429) q[0];
sx q[0];
rz(-1.5786853) q[0];
sx q[0];
rz(-0.042021463) q[0];
rz(0.24629397) q[2];
sx q[2];
rz(-1.8660238) q[2];
sx q[2];
rz(-1.1525046) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.4246108) q[1];
sx q[1];
rz(-1.2094133) q[1];
sx q[1];
rz(-2.7965464) q[1];
x q[2];
rz(-1.7348691) q[3];
sx q[3];
rz(-2.0436156) q[3];
sx q[3];
rz(2.7834322) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.1404672) q[2];
sx q[2];
rz(-0.8320063) q[2];
sx q[2];
rz(2.7531085) q[2];
rz(-1.2515986) q[3];
sx q[3];
rz(-1.3585217) q[3];
sx q[3];
rz(-1.0640915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0851704) q[0];
sx q[0];
rz(-2.6189885) q[0];
sx q[0];
rz(1.1473468) q[0];
rz(-1.7583678) q[1];
sx q[1];
rz(-1.5937832) q[1];
sx q[1];
rz(2.9668818) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6679113) q[0];
sx q[0];
rz(-1.5690049) q[0];
sx q[0];
rz(0.67497336) q[0];
rz(-3.1308451) q[2];
sx q[2];
rz(-0.60950298) q[2];
sx q[2];
rz(-0.94858314) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0510919) q[1];
sx q[1];
rz(-1.8729754) q[1];
sx q[1];
rz(1.7789744) q[1];
x q[2];
rz(-0.21214738) q[3];
sx q[3];
rz(-1.5375932) q[3];
sx q[3];
rz(1.3591059) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.2061578) q[2];
sx q[2];
rz(-1.9400699) q[2];
sx q[2];
rz(-0.60714444) q[2];
rz(-2.855865) q[3];
sx q[3];
rz(-0.61107475) q[3];
sx q[3];
rz(0.40209517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5173335) q[0];
sx q[0];
rz(-1.7835971) q[0];
sx q[0];
rz(0.72878033) q[0];
rz(1.1392611) q[1];
sx q[1];
rz(-1.1092721) q[1];
sx q[1];
rz(1.8702102) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2654289) q[0];
sx q[0];
rz(-3.1294397) q[0];
sx q[0];
rz(-1.7677714) q[0];
rz(-pi) q[1];
x q[1];
rz(0.58614308) q[2];
sx q[2];
rz(-1.8019466) q[2];
sx q[2];
rz(-1.4755206) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.1955586) q[1];
sx q[1];
rz(-0.6094774) q[1];
sx q[1];
rz(-0.35802623) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.0470601) q[3];
sx q[3];
rz(-1.4753398) q[3];
sx q[3];
rz(-0.80194471) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3108958) q[2];
sx q[2];
rz(-1.7252012) q[2];
sx q[2];
rz(2.8719416) q[2];
rz(1.0722748) q[3];
sx q[3];
rz(-0.31338537) q[3];
sx q[3];
rz(2.0872033) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.6894182) q[0];
sx q[0];
rz(-1.9502689) q[0];
sx q[0];
rz(2.7574975) q[0];
rz(-1.1098038) q[1];
sx q[1];
rz(-2.2566819) q[1];
sx q[1];
rz(1.5193411) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.268728) q[0];
sx q[0];
rz(-0.11188398) q[0];
sx q[0];
rz(0.12121339) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.39471925) q[2];
sx q[2];
rz(-2.2119154) q[2];
sx q[2];
rz(1.8973998) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1422524) q[1];
sx q[1];
rz(-0.24079642) q[1];
sx q[1];
rz(-1.1204512) q[1];
x q[2];
rz(-2.8232081) q[3];
sx q[3];
rz(-2.6299713) q[3];
sx q[3];
rz(-1.4871396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.9263837) q[2];
sx q[2];
rz(-2.1314202) q[2];
sx q[2];
rz(2.9877648) q[2];
rz(0.93872968) q[3];
sx q[3];
rz(-1.8925083) q[3];
sx q[3];
rz(1.7083098) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7275823) q[0];
sx q[0];
rz(-1.3879956) q[0];
sx q[0];
rz(-0.54186064) q[0];
rz(-1.2205203) q[1];
sx q[1];
rz(-2.4668756) q[1];
sx q[1];
rz(3.0774679) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.13206088) q[0];
sx q[0];
rz(-1.5144217) q[0];
sx q[0];
rz(1.5359519) q[0];
rz(0.86875963) q[2];
sx q[2];
rz(-1.6693309) q[2];
sx q[2];
rz(-1.2159525) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2623973) q[1];
sx q[1];
rz(-0.91934312) q[1];
sx q[1];
rz(-1.5919973) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7087791) q[3];
sx q[3];
rz(-0.45139957) q[3];
sx q[3];
rz(-1.4572382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.3037783) q[2];
sx q[2];
rz(-2.1393675) q[2];
sx q[2];
rz(-0.079806002) q[2];
rz(-2.5507353) q[3];
sx q[3];
rz(-2.8320524) q[3];
sx q[3];
rz(-3.0965366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87126842) q[0];
sx q[0];
rz(-2.7239983) q[0];
sx q[0];
rz(-2.8840892) q[0];
rz(-0.90266699) q[1];
sx q[1];
rz(-0.84044424) q[1];
sx q[1];
rz(0.88817516) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1290263) q[0];
sx q[0];
rz(-1.2106703) q[0];
sx q[0];
rz(2.9100938) q[0];
rz(-pi) q[1];
x q[1];
rz(0.69239098) q[2];
sx q[2];
rz(-2.2071725) q[2];
sx q[2];
rz(2.6264735) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.17016378) q[1];
sx q[1];
rz(-2.2696643) q[1];
sx q[1];
rz(0.65138833) q[1];
x q[2];
rz(-2.2769663) q[3];
sx q[3];
rz(-2.1451604) q[3];
sx q[3];
rz(0.36353096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.4070134) q[2];
sx q[2];
rz(-2.1740156) q[2];
sx q[2];
rz(0.61545294) q[2];
rz(2.9493799) q[3];
sx q[3];
rz(-2.2902316) q[3];
sx q[3];
rz(-1.6063469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.94083448) q[0];
sx q[0];
rz(-1.2949018) q[0];
sx q[0];
rz(-1.9042263) q[0];
rz(1.0279961) q[1];
sx q[1];
rz(-0.68872394) q[1];
sx q[1];
rz(0.56180305) q[1];
rz(3.0801283) q[2];
sx q[2];
rz(-2.2066084) q[2];
sx q[2];
rz(1.6359272) q[2];
rz(0.54154534) q[3];
sx q[3];
rz(-1.8863746) q[3];
sx q[3];
rz(-2.0538068) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
