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
rz(-2.7849164) q[0];
sx q[0];
rz(-1.4465605) q[0];
sx q[0];
rz(0.90176982) q[0];
rz(4.4737368) q[1];
sx q[1];
rz(2.6982215) q[1];
sx q[1];
rz(8.659521) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1674102) q[0];
sx q[0];
rz(-1.2077792) q[0];
sx q[0];
rz(2.1556751) q[0];
rz(-pi) q[1];
rz(-0.24306007) q[2];
sx q[2];
rz(-1.2420734) q[2];
sx q[2];
rz(-2.6459733) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.81232416) q[1];
sx q[1];
rz(-1.5306889) q[1];
sx q[1];
rz(2.8031111) q[1];
rz(-pi) q[2];
rz(1.8497068) q[3];
sx q[3];
rz(-1.4456914) q[3];
sx q[3];
rz(2.5549614) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.29609933) q[2];
sx q[2];
rz(-2.0534434) q[2];
sx q[2];
rz(1.4720526) q[2];
rz(1.4269525) q[3];
sx q[3];
rz(-1.6441556) q[3];
sx q[3];
rz(2.6078687) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.958309) q[0];
sx q[0];
rz(-1.5077718) q[0];
sx q[0];
rz(-2.2928152) q[0];
rz(-0.035004184) q[1];
sx q[1];
rz(-1.1468381) q[1];
sx q[1];
rz(-0.95357198) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7824421) q[0];
sx q[0];
rz(-1.1138565) q[0];
sx q[0];
rz(0.63553973) q[0];
rz(-pi) q[1];
x q[1];
rz(1.1814647) q[2];
sx q[2];
rz(-1.2876533) q[2];
sx q[2];
rz(2.2623537) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9417022) q[1];
sx q[1];
rz(-0.39559707) q[1];
sx q[1];
rz(-1.0426177) q[1];
x q[2];
rz(-0.93718174) q[3];
sx q[3];
rz(-0.92234334) q[3];
sx q[3];
rz(-2.5536553) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.9780875) q[2];
sx q[2];
rz(-1.3601902) q[2];
sx q[2];
rz(-2.2920442) q[2];
rz(0.16544011) q[3];
sx q[3];
rz(-1.4435507) q[3];
sx q[3];
rz(-0.74976841) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7716832) q[0];
sx q[0];
rz(-2.219438) q[0];
sx q[0];
rz(0.40192303) q[0];
rz(0.15332128) q[1];
sx q[1];
rz(-2.9647277) q[1];
sx q[1];
rz(-2.0053999) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3055666) q[0];
sx q[0];
rz(-1.5587806) q[0];
sx q[0];
rz(-1.1712606) q[0];
x q[1];
rz(-2.6499475) q[2];
sx q[2];
rz(-2.6043713) q[2];
sx q[2];
rz(2.0073089) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.0979251) q[1];
sx q[1];
rz(-0.33722028) q[1];
sx q[1];
rz(-0.80101669) q[1];
x q[2];
rz(0.69092423) q[3];
sx q[3];
rz(-2.5657585) q[3];
sx q[3];
rz(1.3802647) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(3.0574657) q[2];
sx q[2];
rz(-0.87856138) q[2];
sx q[2];
rz(-0.055056661) q[2];
rz(1.347524) q[3];
sx q[3];
rz(-1.3301347) q[3];
sx q[3];
rz(-1.0042892) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.36443001) q[0];
sx q[0];
rz(-0.20474064) q[0];
sx q[0];
rz(-1.5675911) q[0];
rz(0.69152999) q[1];
sx q[1];
rz(-2.4722996) q[1];
sx q[1];
rz(2.9561668) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.90339336) q[0];
sx q[0];
rz(-1.1739587) q[0];
sx q[0];
rz(1.3893886) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.144997) q[2];
sx q[2];
rz(-2.0177207) q[2];
sx q[2];
rz(1.1584692) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.1048829) q[1];
sx q[1];
rz(-1.8817025) q[1];
sx q[1];
rz(-2.8757812) q[1];
rz(-0.33230761) q[3];
sx q[3];
rz(-0.57113591) q[3];
sx q[3];
rz(-0.50779282) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.4039679) q[2];
sx q[2];
rz(-1.4289958) q[2];
sx q[2];
rz(0.0052304012) q[2];
rz(0.38749203) q[3];
sx q[3];
rz(-2.6369075) q[3];
sx q[3];
rz(1.3030049) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
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
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.72652793) q[0];
sx q[0];
rz(-2.1580577) q[0];
sx q[0];
rz(-2.5256185) q[0];
rz(1.9527324) q[1];
sx q[1];
rz(-2.523246) q[1];
sx q[1];
rz(2.628285) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.05840551) q[0];
sx q[0];
rz(-1.8523916) q[0];
sx q[0];
rz(0.0076936184) q[0];
rz(-pi) q[1];
rz(-1.2890069) q[2];
sx q[2];
rz(-2.2766487) q[2];
sx q[2];
rz(-2.6633546) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9037902) q[1];
sx q[1];
rz(-2.2040743) q[1];
sx q[1];
rz(-1.6602181) q[1];
rz(-1.1034455) q[3];
sx q[3];
rz(-1.6765521) q[3];
sx q[3];
rz(0.33973628) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0080228) q[2];
sx q[2];
rz(-1.2761389) q[2];
sx q[2];
rz(-2.3991154) q[2];
rz(-1.6978469) q[3];
sx q[3];
rz(-1.4092813) q[3];
sx q[3];
rz(-1.1013364) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83108574) q[0];
sx q[0];
rz(-2.0948912) q[0];
sx q[0];
rz(2.781784) q[0];
rz(0.85044914) q[1];
sx q[1];
rz(-1.9603739) q[1];
sx q[1];
rz(-2.6757619) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1547928) q[0];
sx q[0];
rz(-1.354711) q[0];
sx q[0];
rz(-2.297202) q[0];
rz(-2.3783422) q[2];
sx q[2];
rz(-1.6692729) q[2];
sx q[2];
rz(-1.5710448) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.1759431) q[1];
sx q[1];
rz(-1.6217983) q[1];
sx q[1];
rz(2.9005364) q[1];
rz(-pi) q[2];
rz(-2.3202002) q[3];
sx q[3];
rz(-0.75858772) q[3];
sx q[3];
rz(2.4241991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.3992074) q[2];
sx q[2];
rz(-1.9611605) q[2];
sx q[2];
rz(-2.7772389) q[2];
rz(2.897701) q[3];
sx q[3];
rz(-0.049592169) q[3];
sx q[3];
rz(-1.2172788) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.57109443) q[0];
sx q[0];
rz(-2.1622393) q[0];
sx q[0];
rz(1.97557) q[0];
rz(-2.1265538) q[1];
sx q[1];
rz(-2.6045585) q[1];
sx q[1];
rz(-0.38645116) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5597661) q[0];
sx q[0];
rz(-2.92972) q[0];
sx q[0];
rz(1.7170402) q[0];
rz(-0.3101686) q[2];
sx q[2];
rz(-1.8156689) q[2];
sx q[2];
rz(-2.819811) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4402953) q[1];
sx q[1];
rz(-3.0616425) q[1];
sx q[1];
rz(-1.246212) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1130969) q[3];
sx q[3];
rz(-1.6784378) q[3];
sx q[3];
rz(0.018317761) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.7361136) q[2];
sx q[2];
rz(-1.4961286) q[2];
sx q[2];
rz(2.9138937) q[2];
rz(-2.0685711) q[3];
sx q[3];
rz(-2.3927972) q[3];
sx q[3];
rz(0.29357114) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88184083) q[0];
sx q[0];
rz(-2.2417534) q[0];
sx q[0];
rz(-2.7135799) q[0];
rz(-2.1233066) q[1];
sx q[1];
rz(-2.7052453) q[1];
sx q[1];
rz(0.87103081) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4081602) q[0];
sx q[0];
rz(-1.0116315) q[0];
sx q[0];
rz(0.59691043) q[0];
x q[1];
rz(1.2134654) q[2];
sx q[2];
rz(-1.7646731) q[2];
sx q[2];
rz(2.9952637) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.5586413) q[1];
sx q[1];
rz(-1.4498596) q[1];
sx q[1];
rz(-1.4287097) q[1];
rz(-pi) q[2];
rz(-0.65555592) q[3];
sx q[3];
rz(-0.79596114) q[3];
sx q[3];
rz(-1.1850921) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.13228664) q[2];
sx q[2];
rz(-2.3460178) q[2];
sx q[2];
rz(-1.2786678) q[2];
rz(-2.5631185) q[3];
sx q[3];
rz(-2.5090802) q[3];
sx q[3];
rz(-1.9875897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4671675) q[0];
sx q[0];
rz(-2.9245057) q[0];
sx q[0];
rz(-2.5417852) q[0];
rz(-1.2990052) q[1];
sx q[1];
rz(-1.5312342) q[1];
sx q[1];
rz(-0.64186796) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0620928) q[0];
sx q[0];
rz(-2.5389517) q[0];
sx q[0];
rz(-0.63061611) q[0];
x q[1];
rz(0.14417837) q[2];
sx q[2];
rz(-0.81071172) q[2];
sx q[2];
rz(-1.7264896) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7781592) q[1];
sx q[1];
rz(-1.5523408) q[1];
sx q[1];
rz(-2.4337971) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5893552) q[3];
sx q[3];
rz(-0.83697666) q[3];
sx q[3];
rz(-2.6058634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.32144) q[2];
sx q[2];
rz(-1.4250616) q[2];
sx q[2];
rz(-2.5313306) q[2];
rz(-2.6214456) q[3];
sx q[3];
rz(-2.0870233) q[3];
sx q[3];
rz(0.49351969) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34058061) q[0];
sx q[0];
rz(-1.8980674) q[0];
sx q[0];
rz(-0.93801671) q[0];
rz(-2.7998789) q[1];
sx q[1];
rz(-1.382117) q[1];
sx q[1];
rz(2.9843073) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.92625916) q[0];
sx q[0];
rz(-2.4865827) q[0];
sx q[0];
rz(2.8192855) q[0];
x q[1];
rz(-1.4862655) q[2];
sx q[2];
rz(-0.8145552) q[2];
sx q[2];
rz(1.5703805) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.966519) q[1];
sx q[1];
rz(-2.0957114) q[1];
sx q[1];
rz(-1.9210084) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1620164) q[3];
sx q[3];
rz(-1.3364387) q[3];
sx q[3];
rz(-2.5046668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.21891317) q[2];
sx q[2];
rz(-2.7965386) q[2];
sx q[2];
rz(1.867713) q[2];
rz(3.1359361) q[3];
sx q[3];
rz(-1.7001245) q[3];
sx q[3];
rz(0.98102942) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(1.2904749) q[0];
sx q[0];
rz(-1.1707476) q[0];
sx q[0];
rz(0.049402417) q[0];
rz(-2.6750917) q[1];
sx q[1];
rz(-1.3460881) q[1];
sx q[1];
rz(0.64429611) q[1];
rz(-0.90032719) q[2];
sx q[2];
rz(-1.1569958) q[2];
sx q[2];
rz(-0.3668084) q[2];
rz(2.7098223) q[3];
sx q[3];
rz(-1.4960066) q[3];
sx q[3];
rz(-3.0638051) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
