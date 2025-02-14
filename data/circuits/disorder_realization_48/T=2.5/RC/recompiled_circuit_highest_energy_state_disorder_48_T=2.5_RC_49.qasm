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
rz(-0.35204044) q[0];
sx q[0];
rz(-0.8249324) q[0];
sx q[0];
rz(-2.6068249) q[0];
rz(-2.1575902) q[1];
sx q[1];
rz(-0.50552955) q[1];
sx q[1];
rz(-1.2619789) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5231294) q[0];
sx q[0];
rz(-2.1697756) q[0];
sx q[0];
rz(1.8904786) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6017351) q[2];
sx q[2];
rz(-1.015402) q[2];
sx q[2];
rz(0.29796165) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.3784005) q[1];
sx q[1];
rz(-2.1175368) q[1];
sx q[1];
rz(-1.631447) q[1];
rz(-2.6747236) q[3];
sx q[3];
rz(-2.7109466) q[3];
sx q[3];
rz(-0.96063411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.24903211) q[2];
sx q[2];
rz(-0.51237115) q[2];
sx q[2];
rz(-1.4830291) q[2];
rz(-0.28462166) q[3];
sx q[3];
rz(-0.62323815) q[3];
sx q[3];
rz(-1.0641789) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2738709) q[0];
sx q[0];
rz(-2.4148648) q[0];
sx q[0];
rz(-3.100585) q[0];
rz(-1.2076591) q[1];
sx q[1];
rz(-2.9137847) q[1];
sx q[1];
rz(-1.6809195) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6448633) q[0];
sx q[0];
rz(-1.71642) q[0];
sx q[0];
rz(-1.54099) q[0];
rz(2.3744205) q[2];
sx q[2];
rz(-1.1572654) q[2];
sx q[2];
rz(-2.3656379) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.25819689) q[1];
sx q[1];
rz(-0.25942311) q[1];
sx q[1];
rz(1.4197465) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3364661) q[3];
sx q[3];
rz(-2.4357492) q[3];
sx q[3];
rz(2.2587551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-0.66551208) q[2];
sx q[2];
rz(-1.6802695) q[2];
sx q[2];
rz(1.3903138) q[2];
rz(-0.0811854) q[3];
sx q[3];
rz(-2.472671) q[3];
sx q[3];
rz(-1.4359052) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7471033) q[0];
sx q[0];
rz(-3.1293226) q[0];
sx q[0];
rz(-2.5315206) q[0];
rz(1.8286145) q[1];
sx q[1];
rz(-0.6956296) q[1];
sx q[1];
rz(-2.5909766) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2064821) q[0];
sx q[0];
rz(-1.2121823) q[0];
sx q[0];
rz(-2.4367843) q[0];
rz(1.2989013) q[2];
sx q[2];
rz(-1.2592794) q[2];
sx q[2];
rz(-2.2443983) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.1296583) q[1];
sx q[1];
rz(-1.7590344) q[1];
sx q[1];
rz(-1.2663575) q[1];
rz(2.9987644) q[3];
sx q[3];
rz(-0.83072829) q[3];
sx q[3];
rz(-2.9252441) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.47784558) q[2];
sx q[2];
rz(-0.17243324) q[2];
sx q[2];
rz(1.0663859) q[2];
rz(-1.1586698) q[3];
sx q[3];
rz(-1.7585157) q[3];
sx q[3];
rz(-3.0995479) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8313507) q[0];
sx q[0];
rz(-0.61315918) q[0];
sx q[0];
rz(-2.2245275) q[0];
rz(-2.4861368) q[1];
sx q[1];
rz(-2.2323445) q[1];
sx q[1];
rz(-0.38958946) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8389908) q[0];
sx q[0];
rz(-0.20379681) q[0];
sx q[0];
rz(1.1558644) q[0];
rz(2.2245313) q[2];
sx q[2];
rz(-0.79981632) q[2];
sx q[2];
rz(0.66674495) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.37871088) q[1];
sx q[1];
rz(-2.7392651) q[1];
sx q[1];
rz(1.4961924) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.54549952) q[3];
sx q[3];
rz(-1.0441458) q[3];
sx q[3];
rz(-0.57034501) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.83170825) q[2];
sx q[2];
rz(-2.1990621) q[2];
sx q[2];
rz(-1.3523098) q[2];
rz(-0.74812198) q[3];
sx q[3];
rz(-1.3040521) q[3];
sx q[3];
rz(-1.6149394) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8378976) q[0];
sx q[0];
rz(-0.56981531) q[0];
sx q[0];
rz(-1.3414398) q[0];
rz(0.98126423) q[1];
sx q[1];
rz(-1.2813247) q[1];
sx q[1];
rz(2.431638) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7959952) q[0];
sx q[0];
rz(-1.5899204) q[0];
sx q[0];
rz(0.32581331) q[0];
x q[1];
rz(2.1469028) q[2];
sx q[2];
rz(-1.369595) q[2];
sx q[2];
rz(-0.76801571) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(0.040474214) q[1];
sx q[1];
rz(-1.3583364) q[1];
sx q[1];
rz(-0.17417769) q[1];
x q[2];
rz(-3.0026765) q[3];
sx q[3];
rz(-0.60231042) q[3];
sx q[3];
rz(1.3773144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.95825163) q[2];
sx q[2];
rz(-2.3812582) q[2];
sx q[2];
rz(2.23488) q[2];
rz(2.9592311) q[3];
sx q[3];
rz(-1.4401108) q[3];
sx q[3];
rz(0.85842925) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7957423) q[0];
sx q[0];
rz(-2.1868571) q[0];
sx q[0];
rz(-0.20172754) q[0];
rz(1.7851104) q[1];
sx q[1];
rz(-2.4701665) q[1];
sx q[1];
rz(1.9904402) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7089823) q[0];
sx q[0];
rz(-1.1469736) q[0];
sx q[0];
rz(2.9124385) q[0];
x q[1];
rz(-2.5207852) q[2];
sx q[2];
rz(-0.70638958) q[2];
sx q[2];
rz(2.8548129) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.21143895) q[1];
sx q[1];
rz(-1.0027998) q[1];
sx q[1];
rz(-1.0878953) q[1];
rz(-1.3428976) q[3];
sx q[3];
rz(-2.0703452) q[3];
sx q[3];
rz(1.7021021) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.98728937) q[2];
sx q[2];
rz(-2.5550948) q[2];
sx q[2];
rz(-2.7941373) q[2];
rz(2.3197428) q[3];
sx q[3];
rz(-1.3653267) q[3];
sx q[3];
rz(1.52389) q[3];
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
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8007941) q[0];
sx q[0];
rz(-0.64685416) q[0];
sx q[0];
rz(1.1232173) q[0];
rz(2.7128291) q[1];
sx q[1];
rz(-2.2518497) q[1];
sx q[1];
rz(-0.14895983) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4892985) q[0];
sx q[0];
rz(-2.6241488) q[0];
sx q[0];
rz(1.6156455) q[0];
rz(-pi) q[1];
rz(-0.2960213) q[2];
sx q[2];
rz(-1.764035) q[2];
sx q[2];
rz(-3.0770965) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.51587869) q[1];
sx q[1];
rz(-0.82059723) q[1];
sx q[1];
rz(-1.3681812) q[1];
rz(-pi) q[2];
rz(1.4895053) q[3];
sx q[3];
rz(-0.40697655) q[3];
sx q[3];
rz(-1.5291884) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.1296156) q[2];
sx q[2];
rz(-0.35795438) q[2];
sx q[2];
rz(-2.820106) q[2];
rz(1.1164411) q[3];
sx q[3];
rz(-0.8588841) q[3];
sx q[3];
rz(-1.4139676) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1534934) q[0];
sx q[0];
rz(-1.3445925) q[0];
sx q[0];
rz(3.1186812) q[0];
rz(1.5494391) q[1];
sx q[1];
rz(-1.0433082) q[1];
sx q[1];
rz(1.2124088) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7717944) q[0];
sx q[0];
rz(-1.8317458) q[0];
sx q[0];
rz(2.1290995) q[0];
x q[1];
rz(-1.9697857) q[2];
sx q[2];
rz(-2.3697457) q[2];
sx q[2];
rz(-3.0400624) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3500994) q[1];
sx q[1];
rz(-1.1731358) q[1];
sx q[1];
rz(-1.2418163) q[1];
rz(-3.0453835) q[3];
sx q[3];
rz(-0.85959496) q[3];
sx q[3];
rz(2.8195153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.2108078) q[2];
sx q[2];
rz(-0.19889861) q[2];
sx q[2];
rz(-0.99736324) q[2];
rz(-1.8831683) q[3];
sx q[3];
rz(-1.1910028) q[3];
sx q[3];
rz(1.2303801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.1579943) q[0];
sx q[0];
rz(-2.4857434) q[0];
sx q[0];
rz(-2.6672145) q[0];
rz(-0.55167088) q[1];
sx q[1];
rz(-2.3730979) q[1];
sx q[1];
rz(2.4840568) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.093134) q[0];
sx q[0];
rz(-3.0952929) q[0];
sx q[0];
rz(-0.53673021) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.9507061) q[2];
sx q[2];
rz(-2.7082241) q[2];
sx q[2];
rz(2.4182408) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.9332757) q[1];
sx q[1];
rz(-0.36872702) q[1];
sx q[1];
rz(1.5255724) q[1];
x q[2];
rz(-0.10513427) q[3];
sx q[3];
rz(-1.7145559) q[3];
sx q[3];
rz(-1.2105699) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-0.23003301) q[2];
sx q[2];
rz(-2.7342789) q[2];
sx q[2];
rz(-1.3113021) q[2];
rz(-2.4274872) q[3];
sx q[3];
rz(-1.8592535) q[3];
sx q[3];
rz(0.56984058) q[3];
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
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.0062200935) q[0];
sx q[0];
rz(-0.94921175) q[0];
sx q[0];
rz(-0.60829341) q[0];
rz(2.3237806) q[1];
sx q[1];
rz(-1.9655971) q[1];
sx q[1];
rz(-2.3086595) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7893164) q[0];
sx q[0];
rz(-2.7171405) q[0];
sx q[0];
rz(0.12571521) q[0];
rz(-pi) q[1];
rz(-1.3141749) q[2];
sx q[2];
rz(-1.2663116) q[2];
sx q[2];
rz(-1.9217132) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.3994308) q[1];
sx q[1];
rz(-0.99282904) q[1];
sx q[1];
rz(-1.4121684) q[1];
rz(-pi) q[2];
rz(0.59304955) q[3];
sx q[3];
rz(-1.5967973) q[3];
sx q[3];
rz(3.0700695) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.1828764) q[2];
sx q[2];
rz(-2.7887838) q[2];
sx q[2];
rz(2.555441) q[2];
rz(2.2385249) q[3];
sx q[3];
rz(-0.62246263) q[3];
sx q[3];
rz(0.085845145) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9027973) q[0];
sx q[0];
rz(-2.0113404) q[0];
sx q[0];
rz(2.197862) q[0];
rz(1.0974274) q[1];
sx q[1];
rz(-0.25249093) q[1];
sx q[1];
rz(2.9846334) q[1];
rz(3.0022754) q[2];
sx q[2];
rz(-1.2323772) q[2];
sx q[2];
rz(2.7215794) q[2];
rz(1.0020574) q[3];
sx q[3];
rz(-2.3593223) q[3];
sx q[3];
rz(0.13025688) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
