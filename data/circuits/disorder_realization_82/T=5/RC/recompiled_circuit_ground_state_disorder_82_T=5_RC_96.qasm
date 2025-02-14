OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.6717186) q[0];
sx q[0];
rz(-1.1550386) q[0];
sx q[0];
rz(1.4299097) q[0];
rz(0.45541304) q[1];
sx q[1];
rz(5.6918511) q[1];
sx q[1];
rz(11.439352) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0394063) q[0];
sx q[0];
rz(-1.8012905) q[0];
sx q[0];
rz(0.011876975) q[0];
rz(-0.65140407) q[2];
sx q[2];
rz(-0.42170516) q[2];
sx q[2];
rz(-2.5159474) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.0309567) q[1];
sx q[1];
rz(-2.1345761) q[1];
sx q[1];
rz(0.33302506) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.5529049) q[3];
sx q[3];
rz(-1.761529) q[3];
sx q[3];
rz(-0.73303661) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0226125) q[2];
sx q[2];
rz(-2.4369414) q[2];
sx q[2];
rz(1.4878147) q[2];
rz(-2.7704499) q[3];
sx q[3];
rz(-2.0149588) q[3];
sx q[3];
rz(2.3625372) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30598518) q[0];
sx q[0];
rz(-1.7636517) q[0];
sx q[0];
rz(0.66407472) q[0];
rz(2.5750419) q[1];
sx q[1];
rz(-1.5357176) q[1];
sx q[1];
rz(-2.9552592) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.60255) q[0];
sx q[0];
rz(-1.3485753) q[0];
sx q[0];
rz(-1.0265539) q[0];
rz(-pi) q[1];
rz(-0.9695942) q[2];
sx q[2];
rz(-2.6896713) q[2];
sx q[2];
rz(2.4659803) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.1809606) q[1];
sx q[1];
rz(-0.69968984) q[1];
sx q[1];
rz(-1.3153428) q[1];
rz(-1.6870895) q[3];
sx q[3];
rz(-2.6175559) q[3];
sx q[3];
rz(2.2268023) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8191007) q[2];
sx q[2];
rz(-1.2995316) q[2];
sx q[2];
rz(1.5619649) q[2];
rz(1.9942079) q[3];
sx q[3];
rz(-0.29427823) q[3];
sx q[3];
rz(1.7779721) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0550302) q[0];
sx q[0];
rz(-1.4921621) q[0];
sx q[0];
rz(0.30558875) q[0];
rz(1.8606404) q[1];
sx q[1];
rz(-2.8260904) q[1];
sx q[1];
rz(1.618128) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9394835) q[0];
sx q[0];
rz(-1.8022707) q[0];
sx q[0];
rz(2.7702296) q[0];
rz(-pi) q[1];
rz(-1.5554713) q[2];
sx q[2];
rz(-1.5256679) q[2];
sx q[2];
rz(-0.082782291) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-3.1393006) q[1];
sx q[1];
rz(-1.2131547) q[1];
sx q[1];
rz(-2.4071818) q[1];
rz(-pi) q[2];
rz(0.95919249) q[3];
sx q[3];
rz(-1.3750031) q[3];
sx q[3];
rz(1.4899561) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.4898701) q[2];
sx q[2];
rz(-1.2469331) q[2];
sx q[2];
rz(2.9215802) q[2];
rz(-0.079719933) q[3];
sx q[3];
rz(-0.97697512) q[3];
sx q[3];
rz(0.93130934) q[3];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43059573) q[0];
sx q[0];
rz(-1.7544704) q[0];
sx q[0];
rz(-0.0084477607) q[0];
rz(-1.9795817) q[1];
sx q[1];
rz(-1.9879257) q[1];
sx q[1];
rz(2.0268424) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.24766391) q[0];
sx q[0];
rz(-1.8801483) q[0];
sx q[0];
rz(0.19750144) q[0];
x q[1];
rz(-0.56416814) q[2];
sx q[2];
rz(-1.086486) q[2];
sx q[2];
rz(-0.51646215) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.1871918) q[1];
sx q[1];
rz(-1.2174509) q[1];
sx q[1];
rz(0.29275972) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.42041789) q[3];
sx q[3];
rz(-1.9501996) q[3];
sx q[3];
rz(0.25015743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.4343425) q[2];
sx q[2];
rz(-1.0184526) q[2];
sx q[2];
rz(0.84558359) q[2];
rz(-0.25104684) q[3];
sx q[3];
rz(-1.8699402) q[3];
sx q[3];
rz(-0.74560753) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.71395981) q[0];
sx q[0];
rz(-0.41249713) q[0];
sx q[0];
rz(-2.1115671) q[0];
rz(2.5469942) q[1];
sx q[1];
rz(-1.4232891) q[1];
sx q[1];
rz(1.3535708) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6973482) q[0];
sx q[0];
rz(-1.591502) q[0];
sx q[0];
rz(0.022788825) q[0];
rz(-pi) q[1];
rz(-1.7556023) q[2];
sx q[2];
rz(-1.2535742) q[2];
sx q[2];
rz(-0.76390195) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.374267) q[1];
sx q[1];
rz(-2.0637636) q[1];
sx q[1];
rz(2.3460991) q[1];
x q[2];
rz(2.5491675) q[3];
sx q[3];
rz(-2.0090143) q[3];
sx q[3];
rz(-0.93417227) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.9109965) q[2];
sx q[2];
rz(-1.5065008) q[2];
sx q[2];
rz(0.050749151) q[2];
rz(2.0283608) q[3];
sx q[3];
rz(-1.166393) q[3];
sx q[3];
rz(-0.39065233) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-2.0473061) q[0];
sx q[0];
rz(-2.081649) q[0];
sx q[0];
rz(2.768709) q[0];
rz(-1.8376384) q[1];
sx q[1];
rz(-1.2645489) q[1];
sx q[1];
rz(1.0252999) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.404843) q[0];
sx q[0];
rz(-1.4944634) q[0];
sx q[0];
rz(0.45053225) q[0];
rz(1.1187024) q[2];
sx q[2];
rz(-2.7774924) q[2];
sx q[2];
rz(0.61227476) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.655788) q[1];
sx q[1];
rz(-0.88690573) q[1];
sx q[1];
rz(2.8832316) q[1];
rz(1.147109) q[3];
sx q[3];
rz(-1.3653339) q[3];
sx q[3];
rz(-0.68483554) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.23124) q[2];
sx q[2];
rz(-2.9515036) q[2];
sx q[2];
rz(-1.7115889) q[2];
rz(1.5405687) q[3];
sx q[3];
rz(-0.68683306) q[3];
sx q[3];
rz(-1.6704667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.98825276) q[0];
sx q[0];
rz(-1.3000458) q[0];
sx q[0];
rz(0.43149313) q[0];
rz(-1.8776114) q[1];
sx q[1];
rz(-1.5411721) q[1];
sx q[1];
rz(-2.4914609) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7509814) q[0];
sx q[0];
rz(-1.6529875) q[0];
sx q[0];
rz(2.85336) q[0];
rz(-pi) q[1];
x q[1];
rz(0.10429485) q[2];
sx q[2];
rz(-2.2513835) q[2];
sx q[2];
rz(-2.7460737) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7909779) q[1];
sx q[1];
rz(-2.2679866) q[1];
sx q[1];
rz(-2.99111) q[1];
rz(-pi) q[2];
x q[2];
rz(1.5863083) q[3];
sx q[3];
rz(-2.6836694) q[3];
sx q[3];
rz(-1.0481038) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5905137) q[2];
sx q[2];
rz(-1.7503909) q[2];
sx q[2];
rz(-0.004465731) q[2];
rz(0.0095857754) q[3];
sx q[3];
rz(-1.8066581) q[3];
sx q[3];
rz(-0.34346223) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.9850605) q[0];
sx q[0];
rz(-2.8592906) q[0];
sx q[0];
rz(2.9811133) q[0];
rz(1.3849244) q[1];
sx q[1];
rz(-2.3424708) q[1];
sx q[1];
rz(0.30881944) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64381448) q[0];
sx q[0];
rz(-0.22499946) q[0];
sx q[0];
rz(-2.1589222) q[0];
rz(1.1868434) q[2];
sx q[2];
rz(-0.89659474) q[2];
sx q[2];
rz(-2.7701258) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.05986) q[1];
sx q[1];
rz(-0.55559056) q[1];
sx q[1];
rz(-2.3742071) q[1];
x q[2];
rz(1.2587955) q[3];
sx q[3];
rz(-1.2332757) q[3];
sx q[3];
rz(2.5520293) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7979692) q[2];
sx q[2];
rz(-0.60911959) q[2];
sx q[2];
rz(1.5999751) q[2];
rz(0.68495098) q[3];
sx q[3];
rz(-1.5870321) q[3];
sx q[3];
rz(1.5233013) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5110382) q[0];
sx q[0];
rz(-1.4383974) q[0];
sx q[0];
rz(2.5901929) q[0];
rz(2.664387) q[1];
sx q[1];
rz(-2.2255032) q[1];
sx q[1];
rz(0.83470693) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3480417) q[0];
sx q[0];
rz(-1.1116331) q[0];
sx q[0];
rz(2.5870917) q[0];
rz(-pi) q[1];
x q[1];
rz(3.1006591) q[2];
sx q[2];
rz(-2.4403095) q[2];
sx q[2];
rz(2.3778039) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1795802) q[1];
sx q[1];
rz(-2.5287345) q[1];
sx q[1];
rz(1.5108766) q[1];
rz(-pi) q[2];
rz(-2.1797997) q[3];
sx q[3];
rz(-2.0941705) q[3];
sx q[3];
rz(-2.9030622) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(1.4386374) q[2];
sx q[2];
rz(-0.96200395) q[2];
sx q[2];
rz(0.23923242) q[2];
rz(2.8171825) q[3];
sx q[3];
rz(-1.7041465) q[3];
sx q[3];
rz(-1.5391866) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5156373) q[0];
sx q[0];
rz(-3.0037168) q[0];
sx q[0];
rz(-0.71664083) q[0];
rz(-0.58249885) q[1];
sx q[1];
rz(-1.7889675) q[1];
sx q[1];
rz(-2.8315721) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4449883) q[0];
sx q[0];
rz(-1.2221081) q[0];
sx q[0];
rz(-0.83072386) q[0];
rz(-pi) q[1];
rz(0.76002174) q[2];
sx q[2];
rz(-2.5993002) q[2];
sx q[2];
rz(-0.063864683) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.9251524) q[1];
sx q[1];
rz(-2.2693949) q[1];
sx q[1];
rz(0.47355814) q[1];
rz(-pi) q[2];
x q[2];
rz(0.72978348) q[3];
sx q[3];
rz(-2.457629) q[3];
sx q[3];
rz(1.3866977) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.6325355) q[2];
sx q[2];
rz(-2.5639503) q[2];
sx q[2];
rz(-1.4198111) q[2];
rz(1.4451197) q[3];
sx q[3];
rz(-0.71075478) q[3];
sx q[3];
rz(2.3783309) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2390908) q[0];
sx q[0];
rz(-0.9466753) q[0];
sx q[0];
rz(-1.7539903) q[0];
rz(-1.6613962) q[1];
sx q[1];
rz(-1.4085242) q[1];
sx q[1];
rz(-2.9396802) q[1];
rz(0.26915941) q[2];
sx q[2];
rz(-2.8648389) q[2];
sx q[2];
rz(1.7493389) q[2];
rz(-0.46764163) q[3];
sx q[3];
rz(-0.94398879) q[3];
sx q[3];
rz(2.332609) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
