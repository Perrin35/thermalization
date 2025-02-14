OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-1.5366323) q[0];
sx q[0];
rz(-0.16083117) q[0];
sx q[0];
rz(-0.39128006) q[0];
rz(-0.04303509) q[1];
sx q[1];
rz(-1.5207091) q[1];
sx q[1];
rz(2.8032141) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7926377) q[0];
sx q[0];
rz(-0.39039877) q[0];
sx q[0];
rz(-0.81116311) q[0];
rz(-pi) q[1];
rz(1.7853043) q[2];
sx q[2];
rz(-1.6202929) q[2];
sx q[2];
rz(-1.6550198) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.0300301) q[1];
sx q[1];
rz(-1.7610368) q[1];
sx q[1];
rz(-0.89221402) q[1];
rz(-pi) q[2];
rz(2.7250763) q[3];
sx q[3];
rz(-1.6702336) q[3];
sx q[3];
rz(0.23867718) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.7732064) q[2];
sx q[2];
rz(-0.97969222) q[2];
sx q[2];
rz(2.0835908) q[2];
rz(0.10231415) q[3];
sx q[3];
rz(-1.5240069) q[3];
sx q[3];
rz(-1.4003096) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31829396) q[0];
sx q[0];
rz(-0.75495356) q[0];
sx q[0];
rz(0.21389432) q[0];
rz(-1.3338044) q[1];
sx q[1];
rz(-0.50024453) q[1];
sx q[1];
rz(2.852829) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.64015111) q[0];
sx q[0];
rz(-2.1955865) q[0];
sx q[0];
rz(-0.21296091) q[0];
rz(-pi) q[1];
rz(-2.3043121) q[2];
sx q[2];
rz(-2.0990666) q[2];
sx q[2];
rz(-2.3579896) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.97280661) q[1];
sx q[1];
rz(-1.3500577) q[1];
sx q[1];
rz(-3.1112063) q[1];
rz(-pi) q[2];
x q[2];
rz(0.64677644) q[3];
sx q[3];
rz(-0.5277165) q[3];
sx q[3];
rz(-2.3257252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-0.6194965) q[2];
sx q[2];
rz(-0.82021004) q[2];
sx q[2];
rz(1.2843457) q[2];
rz(-1.8340825) q[3];
sx q[3];
rz(-1.3176354) q[3];
sx q[3];
rz(-0.6984843) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4035325) q[0];
sx q[0];
rz(-1.0665251) q[0];
sx q[0];
rz(2.2962978) q[0];
rz(-0.90557939) q[1];
sx q[1];
rz(-2.5366492) q[1];
sx q[1];
rz(1.9639429) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2098006) q[0];
sx q[0];
rz(-1.5508979) q[0];
sx q[0];
rz(0.09505247) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.1198172) q[2];
sx q[2];
rz(-2.0284101) q[2];
sx q[2];
rz(2.3624376) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.48124734) q[1];
sx q[1];
rz(-1.3637017) q[1];
sx q[1];
rz(-2.5607263) q[1];
rz(-pi) q[2];
rz(0.92567851) q[3];
sx q[3];
rz(-1.698602) q[3];
sx q[3];
rz(2.5296581) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.5251069) q[2];
sx q[2];
rz(-0.944204) q[2];
sx q[2];
rz(-3.0109829) q[2];
rz(-1.7836001) q[3];
sx q[3];
rz(-1.0087174) q[3];
sx q[3];
rz(1.2880026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.65508771) q[0];
sx q[0];
rz(-2.8772652) q[0];
sx q[0];
rz(2.7274912) q[0];
rz(0.86530238) q[1];
sx q[1];
rz(-1.2369786) q[1];
sx q[1];
rz(0.6032595) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27644409) q[0];
sx q[0];
rz(-1.8521327) q[0];
sx q[0];
rz(-1.9758609) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6418936) q[2];
sx q[2];
rz(-1.9219766) q[2];
sx q[2];
rz(-0.49654135) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.328426) q[1];
sx q[1];
rz(-0.35086497) q[1];
sx q[1];
rz(-3.0163832) q[1];
rz(-pi) q[2];
x q[2];
rz(1.8835265) q[3];
sx q[3];
rz(-1.4566753) q[3];
sx q[3];
rz(-0.097870083) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.52411756) q[2];
sx q[2];
rz(-1.5323428) q[2];
sx q[2];
rz(2.4821846) q[2];
rz(-3.009033) q[3];
sx q[3];
rz(-2.5949251) q[3];
sx q[3];
rz(-2.4373655) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34683126) q[0];
sx q[0];
rz(-2.1857388) q[0];
sx q[0];
rz(2.8723248) q[0];
rz(-1.6412093) q[1];
sx q[1];
rz(-1.2790479) q[1];
sx q[1];
rz(-1.0795116) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74306946) q[0];
sx q[0];
rz(-2.4691104) q[0];
sx q[0];
rz(0.29015707) q[0];
rz(-pi) q[1];
x q[1];
rz(2.3635885) q[2];
sx q[2];
rz(-2.076882) q[2];
sx q[2];
rz(-2.9377459) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.0289356) q[1];
sx q[1];
rz(-1.9017692) q[1];
sx q[1];
rz(2.079936) q[1];
rz(-pi) q[2];
x q[2];
rz(1.7825384) q[3];
sx q[3];
rz(-1.8820888) q[3];
sx q[3];
rz(2.9632729) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.4214804) q[2];
sx q[2];
rz(-0.58028996) q[2];
sx q[2];
rz(-2.271358) q[2];
rz(-1.3358215) q[3];
sx q[3];
rz(-1.3355052) q[3];
sx q[3];
rz(-0.67572063) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0077591) q[0];
sx q[0];
rz(-3.1133911) q[0];
sx q[0];
rz(-2.5303685) q[0];
rz(-2.4080343) q[1];
sx q[1];
rz(-1.4708142) q[1];
sx q[1];
rz(-2.7582817) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4188273) q[0];
sx q[0];
rz(-1.8486029) q[0];
sx q[0];
rz(-0.72577839) q[0];
rz(-pi) q[1];
rz(0.20140751) q[2];
sx q[2];
rz(-2.3039991) q[2];
sx q[2];
rz(-2.7768486) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.94180471) q[1];
sx q[1];
rz(-2.2082303) q[1];
sx q[1];
rz(0.15934039) q[1];
rz(-pi) q[2];
x q[2];
rz(0.46054764) q[3];
sx q[3];
rz(-0.78118284) q[3];
sx q[3];
rz(2.6278842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.7395301) q[2];
sx q[2];
rz(-0.5548839) q[2];
sx q[2];
rz(-2.128672) q[2];
rz(0.5976451) q[3];
sx q[3];
rz(-1.4089855) q[3];
sx q[3];
rz(0.92596084) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0585854) q[0];
sx q[0];
rz(-0.42917955) q[0];
sx q[0];
rz(-0.035932628) q[0];
rz(1.2220471) q[1];
sx q[1];
rz(-0.67600328) q[1];
sx q[1];
rz(-0.24807182) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1897514) q[0];
sx q[0];
rz(-1.3782256) q[0];
sx q[0];
rz(-1.309881) q[0];
rz(-1.7123958) q[2];
sx q[2];
rz(-1.0842825) q[2];
sx q[2];
rz(1.5202735) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.29734941) q[1];
sx q[1];
rz(-1.2037168) q[1];
sx q[1];
rz(-2.966511) q[1];
rz(-pi) q[2];
rz(-0.024240243) q[3];
sx q[3];
rz(-0.10496891) q[3];
sx q[3];
rz(-1.6258383) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.0849358) q[2];
sx q[2];
rz(-2.0818905) q[2];
sx q[2];
rz(-2.6094931) q[2];
rz(-0.45446864) q[3];
sx q[3];
rz(-1.5458509) q[3];
sx q[3];
rz(1.8475629) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4249975) q[0];
sx q[0];
rz(-1.0497365) q[0];
sx q[0];
rz(-2.9014034) q[0];
rz(2.5364618) q[1];
sx q[1];
rz(-1.3477707) q[1];
sx q[1];
rz(-2.4952707) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28531269) q[0];
sx q[0];
rz(-2.3512771) q[0];
sx q[0];
rz(-1.1423443) q[0];
rz(-pi) q[1];
rz(0.036810151) q[2];
sx q[2];
rz(-1.209785) q[2];
sx q[2];
rz(-1.6474468) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.60257116) q[1];
sx q[1];
rz(-0.50316167) q[1];
sx q[1];
rz(-0.43518592) q[1];
x q[2];
rz(-0.21625285) q[3];
sx q[3];
rz(-0.92666221) q[3];
sx q[3];
rz(-2.3210821) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(1.0295082) q[2];
sx q[2];
rz(-1.4936451) q[2];
sx q[2];
rz(3.0599111) q[2];
rz(0.84780848) q[3];
sx q[3];
rz(-0.40926465) q[3];
sx q[3];
rz(-0.22296396) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.34143701) q[0];
sx q[0];
rz(-0.84469047) q[0];
sx q[0];
rz(2.2499114) q[0];
rz(1.2314697) q[1];
sx q[1];
rz(-0.48858085) q[1];
sx q[1];
rz(-0.60985342) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0783892) q[0];
sx q[0];
rz(-1.1097621) q[0];
sx q[0];
rz(1.0151063) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3568809) q[2];
sx q[2];
rz(-1.5810738) q[2];
sx q[2];
rz(2.1267872) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.9111951) q[1];
sx q[1];
rz(-0.814682) q[1];
sx q[1];
rz(-2.7338408) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6173577) q[3];
sx q[3];
rz(-0.34544975) q[3];
sx q[3];
rz(-0.54333401) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.2764728) q[2];
sx q[2];
rz(-0.96729326) q[2];
sx q[2];
rz(2.8709732) q[2];
rz(2.596415) q[3];
sx q[3];
rz(-1.8137648) q[3];
sx q[3];
rz(-2.8934208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4973064) q[0];
sx q[0];
rz(-1.3077284) q[0];
sx q[0];
rz(0.18648952) q[0];
rz(-1.9322152) q[1];
sx q[1];
rz(-1.3554363) q[1];
sx q[1];
rz(0.019651042) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6015897) q[0];
sx q[0];
rz(-2.181291) q[0];
sx q[0];
rz(-2.6363346) q[0];
rz(-pi) q[1];
x q[1];
rz(0.80060963) q[2];
sx q[2];
rz(-0.36663805) q[2];
sx q[2];
rz(-3.0406102) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.8017822) q[1];
sx q[1];
rz(-2.7556917) q[1];
sx q[1];
rz(-2.7245311) q[1];
rz(-pi) q[2];
x q[2];
rz(0.93312937) q[3];
sx q[3];
rz(-2.2266042) q[3];
sx q[3];
rz(0.84049123) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(0.43202117) q[2];
sx q[2];
rz(-2.186128) q[2];
sx q[2];
rz(1.552399) q[2];
rz(-0.10913695) q[3];
sx q[3];
rz(-1.8692632) q[3];
sx q[3];
rz(-2.8652338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
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
rz(-2.6266262) q[0];
sx q[0];
rz(-1.6635386) q[0];
sx q[0];
rz(-1.2141825) q[0];
rz(-1.7264438) q[1];
sx q[1];
rz(-2.1198004) q[1];
sx q[1];
rz(-1.459495) q[1];
rz(-2.696261) q[2];
sx q[2];
rz(-1.8298605) q[2];
sx q[2];
rz(1.4399672) q[2];
rz(0.23701238) q[3];
sx q[3];
rz(-2.7494299) q[3];
sx q[3];
rz(1.3361479) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
