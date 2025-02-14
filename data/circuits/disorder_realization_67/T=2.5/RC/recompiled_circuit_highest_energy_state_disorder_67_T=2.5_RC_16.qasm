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
rz(1.2713852) q[0];
sx q[0];
rz(-0.013590824) q[0];
sx q[0];
rz(3.115227) q[0];
rz(0.68323505) q[1];
sx q[1];
rz(-1.3807715) q[1];
sx q[1];
rz(0.071579054) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.723322) q[0];
sx q[0];
rz(-1.5958565) q[0];
sx q[0];
rz(-1.5805095) q[0];
rz(-0.52710345) q[2];
sx q[2];
rz(-0.44473916) q[2];
sx q[2];
rz(1.2232194) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4178582) q[1];
sx q[1];
rz(-1.5520387) q[1];
sx q[1];
rz(-0.0049519227) q[1];
x q[2];
rz(0.31622131) q[3];
sx q[3];
rz(-1.2653148) q[3];
sx q[3];
rz(-2.9169634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.2804395) q[2];
sx q[2];
rz(-0.70715487) q[2];
sx q[2];
rz(0.46671483) q[2];
rz(-2.6672582) q[3];
sx q[3];
rz(-3.1204087) q[3];
sx q[3];
rz(3.1202313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3050583) q[0];
sx q[0];
rz(-0.49764043) q[0];
sx q[0];
rz(0.019813892) q[0];
rz(1.5927947) q[1];
sx q[1];
rz(-0.21983799) q[1];
sx q[1];
rz(-1.6682495) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23035717) q[0];
sx q[0];
rz(-0.94382554) q[0];
sx q[0];
rz(0.66642739) q[0];
x q[1];
rz(-0.25924087) q[2];
sx q[2];
rz(-1.8564285) q[2];
sx q[2];
rz(-2.223857) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(0.32343601) q[1];
sx q[1];
rz(-1.4940573) q[1];
sx q[1];
rz(1.7652926) q[1];
rz(-2.3109644) q[3];
sx q[3];
rz(-2.1601956) q[3];
sx q[3];
rz(-2.7089861) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-0.3026498) q[2];
sx q[2];
rz(-0.5984211) q[2];
sx q[2];
rz(1.8518651) q[2];
rz(-1.2368115) q[3];
sx q[3];
rz(-2.808282) q[3];
sx q[3];
rz(-0.70241565) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97238338) q[0];
sx q[0];
rz(-1.9820259) q[0];
sx q[0];
rz(-1.5709391) q[0];
rz(1.6698569) q[1];
sx q[1];
rz(-1.5262628) q[1];
sx q[1];
rz(-2.6872046) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.7189302) q[0];
sx q[0];
rz(-0.61394982) q[0];
sx q[0];
rz(1.3387247) q[0];
rz(-pi) q[1];
x q[1];
rz(0.02645385) q[2];
sx q[2];
rz(-1.4288386) q[2];
sx q[2];
rz(0.96645025) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.4207698) q[1];
sx q[1];
rz(-1.5481564) q[1];
sx q[1];
rz(-1.4525502) q[1];
x q[2];
rz(-0.64735712) q[3];
sx q[3];
rz(-2.1212) q[3];
sx q[3];
rz(0.087322012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4001974) q[2];
sx q[2];
rz(-3.1252842) q[2];
sx q[2];
rz(-0.34747094) q[2];
rz(-0.45626429) q[3];
sx q[3];
rz(-0.01472344) q[3];
sx q[3];
rz(-2.1485476) q[3];
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
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4149813) q[0];
sx q[0];
rz(-1.9009637) q[0];
sx q[0];
rz(-1.4062784) q[0];
rz(2.6982488) q[1];
sx q[1];
rz(-1.025238) q[1];
sx q[1];
rz(1.571507) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3298108) q[0];
sx q[0];
rz(-0.64231163) q[0];
sx q[0];
rz(-0.091600939) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.3894044) q[2];
sx q[2];
rz(-3.0471932) q[2];
sx q[2];
rz(1.1122676) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(0.52165495) q[1];
sx q[1];
rz(-1.8490377) q[1];
sx q[1];
rz(-0.010163608) q[1];
rz(1.1394386) q[3];
sx q[3];
rz(-1.316081) q[3];
sx q[3];
rz(1.4875808) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.37322474) q[2];
sx q[2];
rz(-0.39373213) q[2];
sx q[2];
rz(0.079252871) q[2];
rz(1.9521889) q[3];
sx q[3];
rz(-1.3711843) q[3];
sx q[3];
rz(1.6325379) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.70334148) q[0];
sx q[0];
rz(-0.51333135) q[0];
sx q[0];
rz(2.3355423) q[0];
rz(-2.3015859) q[1];
sx q[1];
rz(-3.1286616) q[1];
sx q[1];
rz(-2.3654225) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38473338) q[0];
sx q[0];
rz(-2.1549468) q[0];
sx q[0];
rz(-1.7279051) q[0];
rz(-pi) q[1];
rz(0.010439053) q[2];
sx q[2];
rz(-1.5722646) q[2];
sx q[2];
rz(1.8515585) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.8373903) q[1];
sx q[1];
rz(-1.5588964) q[1];
sx q[1];
rz(-0.13545998) q[1];
rz(1.4587298) q[3];
sx q[3];
rz(-1.8288426) q[3];
sx q[3];
rz(-2.8530981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(3.1157896) q[2];
sx q[2];
rz(-1.5805406) q[2];
sx q[2];
rz(2.4157794) q[2];
rz(-0.18802655) q[3];
sx q[3];
rz(-3.0811716) q[3];
sx q[3];
rz(0.70992011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1793154) q[0];
sx q[0];
rz(-2.582452) q[0];
sx q[0];
rz(2.6002) q[0];
rz(-2.9560282) q[1];
sx q[1];
rz(-1.5488397) q[1];
sx q[1];
rz(3.013179) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5914766) q[0];
sx q[0];
rz(-1.8719561) q[0];
sx q[0];
rz(-2.7635283) q[0];
x q[1];
rz(1.5319848) q[2];
sx q[2];
rz(-0.12141849) q[2];
sx q[2];
rz(-1.613429) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.6278319) q[1];
sx q[1];
rz(-2.7134088) q[1];
sx q[1];
rz(0.19780383) q[1];
x q[2];
rz(-0.010574118) q[3];
sx q[3];
rz(-1.3823783) q[3];
sx q[3];
rz(2.1343975) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7554756) q[2];
sx q[2];
rz(-0.057689276) q[2];
sx q[2];
rz(0.84121394) q[2];
rz(2.9403213) q[3];
sx q[3];
rz(-1.6095251) q[3];
sx q[3];
rz(-0.25963983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5801308) q[0];
sx q[0];
rz(-2.353297) q[0];
sx q[0];
rz(1.5808251) q[0];
rz(0.67281094) q[1];
sx q[1];
rz(-1.7377868) q[1];
sx q[1];
rz(0.028884551) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0085379) q[0];
sx q[0];
rz(-1.8540181) q[0];
sx q[0];
rz(1.3134366) q[0];
x q[1];
rz(-0.37337005) q[2];
sx q[2];
rz(-2.6334116) q[2];
sx q[2];
rz(2.1477107) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.94651949) q[1];
sx q[1];
rz(-1.1596438) q[1];
sx q[1];
rz(2.8821936) q[1];
x q[2];
rz(-2.6050909) q[3];
sx q[3];
rz(-0.85404761) q[3];
sx q[3];
rz(-1.052618) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.2334571) q[2];
sx q[2];
rz(-2.5448749) q[2];
sx q[2];
rz(2.2133568) q[2];
rz(2.8440031) q[3];
sx q[3];
rz(-2.987515) q[3];
sx q[3];
rz(0.84460622) q[3];
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
rz(-3.0723648) q[0];
sx q[0];
rz(-0.2437676) q[0];
sx q[0];
rz(0.099076554) q[0];
rz(0.98723269) q[1];
sx q[1];
rz(-1.3039373) q[1];
sx q[1];
rz(-0.5388906) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19463704) q[0];
sx q[0];
rz(-1.6536923) q[0];
sx q[0];
rz(-3.1350053) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6201934) q[2];
sx q[2];
rz(-1.5461174) q[2];
sx q[2];
rz(-0.53804735) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.1325593) q[1];
sx q[1];
rz(-1.0485639) q[1];
sx q[1];
rz(1.9886677) q[1];
rz(-0.022607443) q[3];
sx q[3];
rz(-1.5033098) q[3];
sx q[3];
rz(2.8565503) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.99523669) q[2];
sx q[2];
rz(-0.015492798) q[2];
sx q[2];
rz(-2.8323979) q[2];
rz(2.501798) q[3];
sx q[3];
rz(-0.00034172405) q[3];
sx q[3];
rz(3.0777847) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.30534202) q[0];
sx q[0];
rz(-0.59665614) q[0];
sx q[0];
rz(-0.054656595) q[0];
rz(1.1326185) q[1];
sx q[1];
rz(-1.9839958) q[1];
sx q[1];
rz(1.2861015) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.085145935) q[0];
sx q[0];
rz(-2.9582199) q[0];
sx q[0];
rz(1.8118748) q[0];
x q[1];
rz(3.1007721) q[2];
sx q[2];
rz(-1.5287182) q[2];
sx q[2];
rz(0.034860858) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.1117592) q[1];
sx q[1];
rz(-2.8984424) q[1];
sx q[1];
rz(0.38317005) q[1];
rz(-pi) q[2];
rz(-1.4350227) q[3];
sx q[3];
rz(-2.6505436) q[3];
sx q[3];
rz(0.74573475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.9201811) q[2];
sx q[2];
rz(-0.56500089) q[2];
sx q[2];
rz(1.1037702) q[2];
rz(-1.6086027) q[3];
sx q[3];
rz(-0.04403232) q[3];
sx q[3];
rz(-0.65582961) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.00103818) q[0];
sx q[0];
rz(-0.1796722) q[0];
sx q[0];
rz(-3.1371064) q[0];
rz(1.5525612) q[1];
sx q[1];
rz(-1.4493425) q[1];
sx q[1];
rz(0.057131279) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73308459) q[0];
sx q[0];
rz(-1.391469) q[0];
sx q[0];
rz(-3.0441557) q[0];
x q[1];
rz(-0.22259076) q[2];
sx q[2];
rz(-1.711297) q[2];
sx q[2];
rz(-2.8369571) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.88567121) q[1];
sx q[1];
rz(-2.2388865) q[1];
sx q[1];
rz(-2.9873965) q[1];
rz(-pi) q[2];
x q[2];
rz(0.1520429) q[3];
sx q[3];
rz(-0.83629823) q[3];
sx q[3];
rz(-1.9805857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.750018) q[2];
sx q[2];
rz(-0.027316814) q[2];
sx q[2];
rz(0.86891437) q[2];
rz(-0.9817552) q[3];
sx q[3];
rz(-3.1120286) q[3];
sx q[3];
rz(0.50299197) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0960196) q[0];
sx q[0];
rz(-1.6572784) q[0];
sx q[0];
rz(-1.4838765) q[0];
rz(2.6847196) q[1];
sx q[1];
rz(-0.15468205) q[1];
sx q[1];
rz(-0.044943132) q[1];
rz(-0.1931242) q[2];
sx q[2];
rz(-2.3681691) q[2];
sx q[2];
rz(-2.9409627) q[2];
rz(1.7140688) q[3];
sx q[3];
rz(-2.2289056) q[3];
sx q[3];
rz(-3.0994305) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
