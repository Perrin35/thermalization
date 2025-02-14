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
rz(0.96639291) q[0];
sx q[0];
rz(-1.7557431) q[0];
sx q[0];
rz(2.5461499) q[0];
rz(-1.8499941) q[1];
sx q[1];
rz(-0.59102494) q[1];
sx q[1];
rz(0.12330595) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.3759252) q[0];
sx q[0];
rz(-2.1095089) q[0];
sx q[0];
rz(-0.49077351) q[0];
rz(-pi) q[1];
x q[1];
rz(0.79084556) q[2];
sx q[2];
rz(-0.95480761) q[2];
sx q[2];
rz(-3.0012584) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.6001624) q[1];
sx q[1];
rz(-1.7847536) q[1];
sx q[1];
rz(-3.0453277) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7756157) q[3];
sx q[3];
rz(-1.4282188) q[3];
sx q[3];
rz(-2.3516084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.3143602) q[2];
sx q[2];
rz(-2.0388956) q[2];
sx q[2];
rz(3.0296791) q[2];
rz(1.7498451) q[3];
sx q[3];
rz(-1.9839957) q[3];
sx q[3];
rz(-0.30645034) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7258485) q[0];
sx q[0];
rz(-2.2968676) q[0];
sx q[0];
rz(-2.9916905) q[0];
rz(2.8556178) q[1];
sx q[1];
rz(-1.3949225) q[1];
sx q[1];
rz(2.012595) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1492831) q[0];
sx q[0];
rz(-2.8516927) q[0];
sx q[0];
rz(-1.4647746) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.01458077) q[2];
sx q[2];
rz(-1.5111856) q[2];
sx q[2];
rz(-2.4042442) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.2098238) q[1];
sx q[1];
rz(-1.2126414) q[1];
sx q[1];
rz(1.7478554) q[1];
x q[2];
rz(-2.1993447) q[3];
sx q[3];
rz(-2.0671004) q[3];
sx q[3];
rz(-1.7402349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.67843208) q[2];
sx q[2];
rz(-2.2191935) q[2];
sx q[2];
rz(1.9909667) q[2];
rz(1.1848508) q[3];
sx q[3];
rz(-1.4169644) q[3];
sx q[3];
rz(3.1209893) q[3];
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
rz(-pi) q[0];
sx q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6612369) q[0];
sx q[0];
rz(-2.895597) q[0];
sx q[0];
rz(2.7599957) q[0];
rz(1.2941788) q[1];
sx q[1];
rz(-1.1062063) q[1];
sx q[1];
rz(0.19827422) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0509153) q[0];
sx q[0];
rz(-2.8534575) q[0];
sx q[0];
rz(-2.227562) q[0];
rz(-pi) q[1];
rz(-2.2402693) q[2];
sx q[2];
rz(-1.5076037) q[2];
sx q[2];
rz(-2.2128999) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.72858649) q[1];
sx q[1];
rz(-2.8236514) q[1];
sx q[1];
rz(2.4012662) q[1];
rz(-pi) q[2];
x q[2];
rz(1.4573077) q[3];
sx q[3];
rz(-0.17593613) q[3];
sx q[3];
rz(-0.709155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.3759489) q[2];
sx q[2];
rz(-2.9889034) q[2];
sx q[2];
rz(0.51885968) q[2];
rz(-1.8910003) q[3];
sx q[3];
rz(-1.0873245) q[3];
sx q[3];
rz(-0.58108228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4984703) q[0];
sx q[0];
rz(-2.9415218) q[0];
sx q[0];
rz(2.6923687) q[0];
rz(-1.6953702) q[1];
sx q[1];
rz(-0.64165533) q[1];
sx q[1];
rz(0.62072388) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.45263824) q[0];
sx q[0];
rz(-1.7151378) q[0];
sx q[0];
rz(2.4792433) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.63346699) q[2];
sx q[2];
rz(-1.4711498) q[2];
sx q[2];
rz(-1.5758621) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.6606969) q[1];
sx q[1];
rz(-0.45240739) q[1];
sx q[1];
rz(-0.47641944) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4820547) q[3];
sx q[3];
rz(-2.6013298) q[3];
sx q[3];
rz(2.1335354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.1064523) q[2];
sx q[2];
rz(-1.8702714) q[2];
sx q[2];
rz(-0.16109666) q[2];
rz(2.7437239) q[3];
sx q[3];
rz(-1.4784003) q[3];
sx q[3];
rz(0.14959344) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.38254) q[0];
sx q[0];
rz(-1.362514) q[0];
sx q[0];
rz(0.25704849) q[0];
rz(2.2611179) q[1];
sx q[1];
rz(-0.6404666) q[1];
sx q[1];
rz(1.2535198) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3220164) q[0];
sx q[0];
rz(-1.5347693) q[0];
sx q[0];
rz(-3.1234804) q[0];
rz(2.9078616) q[2];
sx q[2];
rz(-2.3178551) q[2];
sx q[2];
rz(1.9875789) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.94379986) q[1];
sx q[1];
rz(-1.5256923) q[1];
sx q[1];
rz(2.3438441) q[1];
rz(1.5725617) q[3];
sx q[3];
rz(-1.2236508) q[3];
sx q[3];
rz(2.7070759) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.72838655) q[2];
sx q[2];
rz(-1.2089968) q[2];
sx q[2];
rz(-1.8750635) q[2];
rz(2.5201216) q[3];
sx q[3];
rz(-2.5340243) q[3];
sx q[3];
rz(0.5411886) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.40015873) q[0];
sx q[0];
rz(-2.3024237) q[0];
sx q[0];
rz(-0.025064502) q[0];
rz(1.9301682) q[1];
sx q[1];
rz(-1.8823267) q[1];
sx q[1];
rz(-2.8866344) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4716041) q[0];
sx q[0];
rz(-1.5657664) q[0];
sx q[0];
rz(0.049335376) q[0];
x q[1];
rz(1.0827052) q[2];
sx q[2];
rz(-1.0399264) q[2];
sx q[2];
rz(-1.7386029) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3759758) q[1];
sx q[1];
rz(-1.0436397) q[1];
sx q[1];
rz(-0.65816452) q[1];
rz(2.1433349) q[3];
sx q[3];
rz(-1.7504331) q[3];
sx q[3];
rz(1.6999753) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.8956464) q[2];
sx q[2];
rz(-1.0422372) q[2];
sx q[2];
rz(-0.45905217) q[2];
rz(1.6445271) q[3];
sx q[3];
rz(-0.11681695) q[3];
sx q[3];
rz(-3.0204401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46886214) q[0];
sx q[0];
rz(-2.6519096) q[0];
sx q[0];
rz(0.086061867) q[0];
rz(-2.22279) q[1];
sx q[1];
rz(-2.0252392) q[1];
sx q[1];
rz(-1.2109717) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.051832599) q[0];
sx q[0];
rz(-2.7564232) q[0];
sx q[0];
rz(-2.5776107) q[0];
rz(-2.7007799) q[2];
sx q[2];
rz(-2.0939504) q[2];
sx q[2];
rz(2.3135106) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.705928) q[1];
sx q[1];
rz(-1.8907232) q[1];
sx q[1];
rz(-2.9271057) q[1];
x q[2];
rz(-2.7769412) q[3];
sx q[3];
rz(-1.9130008) q[3];
sx q[3];
rz(2.3221644) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.81898895) q[2];
sx q[2];
rz(-2.6456867) q[2];
sx q[2];
rz(-0.71713478) q[2];
rz(1.0273733) q[3];
sx q[3];
rz(-1.7337948) q[3];
sx q[3];
rz(0.43924847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1502007) q[0];
sx q[0];
rz(-2.1450277) q[0];
sx q[0];
rz(-3.126934) q[0];
rz(-2.390059) q[1];
sx q[1];
rz(-1.9207759) q[1];
sx q[1];
rz(-1.6709447) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.178407) q[0];
sx q[0];
rz(-1.3826332) q[0];
sx q[0];
rz(1.7117731) q[0];
rz(-pi) q[1];
rz(1.650943) q[2];
sx q[2];
rz(-1.8266457) q[2];
sx q[2];
rz(0.7776153) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4336509) q[1];
sx q[1];
rz(-1.4580863) q[1];
sx q[1];
rz(-1.9757084) q[1];
x q[2];
rz(0.003753547) q[3];
sx q[3];
rz(-1.9158746) q[3];
sx q[3];
rz(-1.0964637) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.8810001) q[2];
sx q[2];
rz(-1.1129817) q[2];
sx q[2];
rz(0.27349681) q[2];
rz(-1.986844) q[3];
sx q[3];
rz(-2.4879849) q[3];
sx q[3];
rz(2.4912513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.038789373) q[0];
sx q[0];
rz(-1.1253072) q[0];
sx q[0];
rz(-1.0754841) q[0];
rz(0.12818809) q[1];
sx q[1];
rz(-2.2905541) q[1];
sx q[1];
rz(-0.34559616) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.95042607) q[0];
sx q[0];
rz(-1.4109525) q[0];
sx q[0];
rz(2.2402083) q[0];
x q[1];
rz(2.0515726) q[2];
sx q[2];
rz(-0.62453237) q[2];
sx q[2];
rz(-2.9337954) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.83452713) q[1];
sx q[1];
rz(-2.5967) q[1];
sx q[1];
rz(-0.12172555) q[1];
x q[2];
rz(2.9008222) q[3];
sx q[3];
rz(-0.61175377) q[3];
sx q[3];
rz(1.4021378) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.0547611) q[2];
sx q[2];
rz(-2.0589477) q[2];
sx q[2];
rz(-1.6205988) q[2];
rz(1.7704376) q[3];
sx q[3];
rz(-2.3833279) q[3];
sx q[3];
rz(0.41485205) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3138251) q[0];
sx q[0];
rz(-0.96243745) q[0];
sx q[0];
rz(-2.9058822) q[0];
rz(0.57304263) q[1];
sx q[1];
rz(-1.0083116) q[1];
sx q[1];
rz(-1.7726353) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0018250759) q[0];
sx q[0];
rz(-2.2968676) q[0];
sx q[0];
rz(2.1609309) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1718352) q[2];
sx q[2];
rz(-2.1548993) q[2];
sx q[2];
rz(0.36512185) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.072242529) q[1];
sx q[1];
rz(-1.4003039) q[1];
sx q[1];
rz(-2.1428698) q[1];
rz(-pi) q[2];
x q[2];
rz(0.36054109) q[3];
sx q[3];
rz(-1.6586132) q[3];
sx q[3];
rz(-2.5537864) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3706751) q[2];
sx q[2];
rz(-1.3266027) q[2];
sx q[2];
rz(0.053000432) q[2];
rz(-0.75183374) q[3];
sx q[3];
rz(-2.1078608) q[3];
sx q[3];
rz(-2.2161765) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5175405) q[0];
sx q[0];
rz(-2.1255827) q[0];
sx q[0];
rz(-0.95638635) q[0];
rz(1.3507631) q[1];
sx q[1];
rz(-0.39719926) q[1];
sx q[1];
rz(-0.15695708) q[1];
rz(-0.69722947) q[2];
sx q[2];
rz(-1.4705428) q[2];
sx q[2];
rz(-1.85208) q[2];
rz(-2.1696321) q[3];
sx q[3];
rz(-2.2836015) q[3];
sx q[3];
rz(-2.9290269) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
