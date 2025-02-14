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
rz(0.056869153) q[0];
sx q[0];
rz(-0.19357227) q[0];
sx q[0];
rz(-0.40785664) q[0];
rz(-0.017539311) q[1];
sx q[1];
rz(4.3932228) q[1];
sx q[1];
rz(7.8235758) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2112132) q[0];
sx q[0];
rz(-1.6193074) q[0];
sx q[0];
rz(-2.0082974) q[0];
rz(1.4354188) q[2];
sx q[2];
rz(-1.1203566) q[2];
sx q[2];
rz(3.0036486) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(1.7471024) q[1];
sx q[1];
rz(-1.1848406) q[1];
sx q[1];
rz(-1.3420022) q[1];
rz(-pi) q[2];
rz(2.3204221) q[3];
sx q[3];
rz(-1.0945012) q[3];
sx q[3];
rz(1.1007635) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.5876329) q[2];
sx q[2];
rz(-2.0822058) q[2];
sx q[2];
rz(-0.21767347) q[2];
rz(-0.31146464) q[3];
sx q[3];
rz(-0.61190999) q[3];
sx q[3];
rz(-1.9757087) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.72164732) q[0];
sx q[0];
rz(-0.2758652) q[0];
sx q[0];
rz(-2.4496147) q[0];
rz(0.75633374) q[1];
sx q[1];
rz(-2.6192009) q[1];
sx q[1];
rz(-1.5501032) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.67705757) q[0];
sx q[0];
rz(-1.4486396) q[0];
sx q[0];
rz(0.71420251) q[0];
rz(1.1728376) q[2];
sx q[2];
rz(-1.6670456) q[2];
sx q[2];
rz(-2.725696) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.8988425) q[1];
sx q[1];
rz(-1.4489343) q[1];
sx q[1];
rz(-0.16074462) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.85568537) q[3];
sx q[3];
rz(-0.69310235) q[3];
sx q[3];
rz(-0.86831304) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.1456566) q[2];
sx q[2];
rz(-2.9684976) q[2];
sx q[2];
rz(0.23925979) q[2];
rz(1.4907106) q[3];
sx q[3];
rz(-1.2942634) q[3];
sx q[3];
rz(1.9132805) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2759129) q[0];
sx q[0];
rz(-0.90621197) q[0];
sx q[0];
rz(-3.1269585) q[0];
rz(-2.2485661) q[1];
sx q[1];
rz(-2.1664186) q[1];
sx q[1];
rz(-0.21569529) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.53447002) q[0];
sx q[0];
rz(-1.5661582) q[0];
sx q[0];
rz(1.5195373) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.2273618) q[2];
sx q[2];
rz(-0.94007713) q[2];
sx q[2];
rz(-1.3892) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.5657422) q[1];
sx q[1];
rz(-1.3389968) q[1];
sx q[1];
rz(-1.4398458) q[1];
rz(2.2641597) q[3];
sx q[3];
rz(-1.0022396) q[3];
sx q[3];
rz(1.3662149) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-1.0594242) q[2];
sx q[2];
rz(-1.1801722) q[2];
sx q[2];
rz(-2.1480985) q[2];
rz(-0.70837402) q[3];
sx q[3];
rz(-1.1185442) q[3];
sx q[3];
rz(-0.12349252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3941536) q[0];
sx q[0];
rz(-2.6469632) q[0];
sx q[0];
rz(-0.69394773) q[0];
rz(-2.8548062) q[1];
sx q[1];
rz(-1.5319805) q[1];
sx q[1];
rz(2.0538816) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74831731) q[0];
sx q[0];
rz(-1.4469742) q[0];
sx q[0];
rz(2.078915) q[0];
rz(-pi) q[1];
x q[1];
rz(2.094628) q[2];
sx q[2];
rz(-2.0598542) q[2];
sx q[2];
rz(-0.2917052) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.4162594) q[1];
sx q[1];
rz(-0.97914588) q[1];
sx q[1];
rz(0.36794956) q[1];
x q[2];
rz(0.71432738) q[3];
sx q[3];
rz(-1.6052402) q[3];
sx q[3];
rz(-1.2726651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9168758) q[2];
sx q[2];
rz(-1.1563053) q[2];
sx q[2];
rz(1.7737596) q[2];
rz(-1.9035089) q[3];
sx q[3];
rz(-1.5690469) q[3];
sx q[3];
rz(2.9253173) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
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
rz(1.5424159) q[0];
sx q[0];
rz(-1.8734064) q[0];
sx q[0];
rz(-2.5237778) q[0];
rz(-1.1082331) q[1];
sx q[1];
rz(-0.30312678) q[1];
sx q[1];
rz(-2.2584426) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9734479) q[0];
sx q[0];
rz(-0.75605481) q[0];
sx q[0];
rz(3.130828) q[0];
x q[1];
rz(-0.43608433) q[2];
sx q[2];
rz(-1.423342) q[2];
sx q[2];
rz(2.6718692) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.0269491) q[1];
sx q[1];
rz(-1.0806688) q[1];
sx q[1];
rz(-3.0672706) q[1];
rz(-pi) q[2];
x q[2];
rz(0.95728504) q[3];
sx q[3];
rz(-2.4924879) q[3];
sx q[3];
rz(1.7136991) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.9250179) q[2];
sx q[2];
rz(-0.25187945) q[2];
sx q[2];
rz(1.7803171) q[2];
rz(-1.2733634) q[3];
sx q[3];
rz(-1.4342156) q[3];
sx q[3];
rz(2.1586965) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-0.044416044) q[0];
sx q[0];
rz(-1.2460848) q[0];
sx q[0];
rz(0.62028766) q[0];
rz(2.7445131) q[1];
sx q[1];
rz(-2.670791) q[1];
sx q[1];
rz(-1.9420067) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.37017779) q[0];
sx q[0];
rz(-1.6770898) q[0];
sx q[0];
rz(-3.1205157) q[0];
rz(-pi) q[1];
rz(1.7825837) q[2];
sx q[2];
rz(-1.7149936) q[2];
sx q[2];
rz(1.5825001) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.6288249) q[1];
sx q[1];
rz(-0.56124748) q[1];
sx q[1];
rz(0.36833771) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.32034822) q[3];
sx q[3];
rz(-0.36133063) q[3];
sx q[3];
rz(0.02905127) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8951796) q[2];
sx q[2];
rz(-1.8975039) q[2];
sx q[2];
rz(-0.67071521) q[2];
rz(-0.29117584) q[3];
sx q[3];
rz(-0.61101919) q[3];
sx q[3];
rz(0.62180716) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.42596844) q[0];
sx q[0];
rz(-1.427303) q[0];
sx q[0];
rz(-0.17844644) q[0];
rz(-2.2715691) q[1];
sx q[1];
rz(-0.88302892) q[1];
sx q[1];
rz(-2.3387486) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.298379) q[0];
sx q[0];
rz(-1.5456539) q[0];
sx q[0];
rz(-2.521994) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.7807671) q[2];
sx q[2];
rz(-0.61491167) q[2];
sx q[2];
rz(-0.71984839) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.1969779) q[1];
sx q[1];
rz(-0.77003819) q[1];
sx q[1];
rz(0.47702392) q[1];
rz(-pi) q[2];
rz(-1.22193) q[3];
sx q[3];
rz(-2.4838944) q[3];
sx q[3];
rz(-1.5869092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.37658438) q[2];
sx q[2];
rz(-1.2976126) q[2];
sx q[2];
rz(-0.89361781) q[2];
rz(-1.5417967) q[3];
sx q[3];
rz(-1.4897646) q[3];
sx q[3];
rz(-2.5347575) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1641418) q[0];
sx q[0];
rz(-0.41005382) q[0];
sx q[0];
rz(0.2151016) q[0];
rz(2.2970301) q[1];
sx q[1];
rz(-1.8212049) q[1];
sx q[1];
rz(0.79427687) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0117019) q[0];
sx q[0];
rz(-0.029057682) q[0];
sx q[0];
rz(-0.41883166) q[0];
rz(-pi) q[1];
rz(-1.2634981) q[2];
sx q[2];
rz(-1.1329044) q[2];
sx q[2];
rz(1.7192993) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.46103288) q[1];
sx q[1];
rz(-2.5941412) q[1];
sx q[1];
rz(1.7458818) q[1];
rz(-2.8135962) q[3];
sx q[3];
rz(-2.3324926) q[3];
sx q[3];
rz(-2.3878333) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(0.27791417) q[2];
sx q[2];
rz(-1.2028376) q[2];
sx q[2];
rz(1.4097144) q[2];
rz(-2.388741) q[3];
sx q[3];
rz(-1.1466305) q[3];
sx q[3];
rz(1.0021771) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.79733217) q[0];
sx q[0];
rz(-0.40818885) q[0];
sx q[0];
rz(0.97287384) q[0];
rz(-1.6196039) q[1];
sx q[1];
rz(-1.7771959) q[1];
sx q[1];
rz(2.5111228) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3909633) q[0];
sx q[0];
rz(-1.7087775) q[0];
sx q[0];
rz(-0.18355455) q[0];
rz(-pi) q[1];
rz(1.5233598) q[2];
sx q[2];
rz(-2.9844797) q[2];
sx q[2];
rz(-0.6873129) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.7122102) q[1];
sx q[1];
rz(-2.3377067) q[1];
sx q[1];
rz(-1.5014807) q[1];
x q[2];
rz(0.34548605) q[3];
sx q[3];
rz(-2.957323) q[3];
sx q[3];
rz(0.94217448) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-3.1260599) q[2];
sx q[2];
rz(-1.9893179) q[2];
sx q[2];
rz(1.3014334) q[2];
rz(1.556501) q[3];
sx q[3];
rz(-2.3157178) q[3];
sx q[3];
rz(-2.7263156) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.62366098) q[0];
sx q[0];
rz(-0.15075891) q[0];
sx q[0];
rz(-0.49322042) q[0];
rz(-1.2552235) q[1];
sx q[1];
rz(-1.8938277) q[1];
sx q[1];
rz(0.36466041) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5887774) q[0];
sx q[0];
rz(-1.7087666) q[0];
sx q[0];
rz(1.9068043) q[0];
rz(-pi) q[1];
rz(2.2232391) q[2];
sx q[2];
rz(-1.9404354) q[2];
sx q[2];
rz(0.015394848) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.11110567) q[1];
sx q[1];
rz(-2.3685622) q[1];
sx q[1];
rz(1.4902601) q[1];
x q[2];
rz(1.4480035) q[3];
sx q[3];
rz(-2.8229575) q[3];
sx q[3];
rz(2.2991869) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.4517335) q[2];
sx q[2];
rz(-1.006459) q[2];
sx q[2];
rz(0.9922007) q[2];
rz(-0.36710468) q[3];
sx q[3];
rz(-1.7016943) q[3];
sx q[3];
rz(-0.67830694) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.622396) q[0];
sx q[0];
rz(-1.4764897) q[0];
sx q[0];
rz(1.3523703) q[0];
rz(-0.92338152) q[1];
sx q[1];
rz(-2.2877749) q[1];
sx q[1];
rz(1.9705082) q[1];
rz(2.1869356) q[2];
sx q[2];
rz(-2.8980394) q[2];
sx q[2];
rz(0.63799636) q[2];
rz(2.4319875) q[3];
sx q[3];
rz(-0.81021877) q[3];
sx q[3];
rz(2.6859663) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
