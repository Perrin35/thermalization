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
rz(-0.77684075) q[0];
sx q[0];
rz(-0.87640327) q[0];
sx q[0];
rz(-0.25282282) q[0];
rz(1.1480992) q[1];
sx q[1];
rz(-2.2263081) q[1];
sx q[1];
rz(-0.80438703) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0641805) q[0];
sx q[0];
rz(-1.8837116) q[0];
sx q[0];
rz(1.1699647) q[0];
rz(-1.3374653) q[2];
sx q[2];
rz(-1.1510282) q[2];
sx q[2];
rz(0.45316089) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.3004669) q[1];
sx q[1];
rz(-1.9411095) q[1];
sx q[1];
rz(-1.200586) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.4553045) q[3];
sx q[3];
rz(-1.9126429) q[3];
sx q[3];
rz(-1.6153743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.048451) q[2];
sx q[2];
rz(-2.1790049) q[2];
sx q[2];
rz(0.87240458) q[2];
rz(0.33556542) q[3];
sx q[3];
rz(-1.7458956) q[3];
sx q[3];
rz(-3.057835) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63808477) q[0];
sx q[0];
rz(-0.26569772) q[0];
sx q[0];
rz(0.22915325) q[0];
rz(-0.19042641) q[1];
sx q[1];
rz(-1.8016022) q[1];
sx q[1];
rz(0.71358877) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.23799831) q[0];
sx q[0];
rz(-2.2269214) q[0];
sx q[0];
rz(2.800992) q[0];
rz(0.53411463) q[2];
sx q[2];
rz(-2.2306109) q[2];
sx q[2];
rz(0.5420064) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.16878) q[1];
sx q[1];
rz(-2.3832364) q[1];
sx q[1];
rz(-2.7060164) q[1];
x q[2];
rz(1.7341033) q[3];
sx q[3];
rz(-2.3288832) q[3];
sx q[3];
rz(0.14369609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.4552292) q[2];
sx q[2];
rz(-2.8296622) q[2];
sx q[2];
rz(-2.6658106) q[2];
rz(-1.0698211) q[3];
sx q[3];
rz(-1.0082303) q[3];
sx q[3];
rz(0.76655918) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0723202) q[0];
sx q[0];
rz(-1.197149) q[0];
sx q[0];
rz(-1.9092165) q[0];
rz(2.6909289) q[1];
sx q[1];
rz(-1.181299) q[1];
sx q[1];
rz(-0.88952363) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33143932) q[0];
sx q[0];
rz(-0.74140775) q[0];
sx q[0];
rz(-2.5196694) q[0];
rz(-pi) q[1];
x q[1];
rz(0.81368229) q[2];
sx q[2];
rz(-2.431834) q[2];
sx q[2];
rz(-2.8738632) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.079021318) q[1];
sx q[1];
rz(-1.3600742) q[1];
sx q[1];
rz(-0.91075588) q[1];
rz(2.7660288) q[3];
sx q[3];
rz(-1.8960388) q[3];
sx q[3];
rz(2.2328245) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.4543317) q[2];
sx q[2];
rz(-0.4250409) q[2];
sx q[2];
rz(-1.8827776) q[2];
rz(1.2339633) q[3];
sx q[3];
rz(-1.1326658) q[3];
sx q[3];
rz(3.1335355) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
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
rz(-2.4312129) q[0];
sx q[0];
rz(-2.2358535) q[0];
sx q[0];
rz(1.4696962) q[0];
rz(-1.1471033) q[1];
sx q[1];
rz(-1.0018145) q[1];
sx q[1];
rz(-2.240644) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1484447) q[0];
sx q[0];
rz(-1.6560418) q[0];
sx q[0];
rz(-2.8277525) q[0];
x q[1];
rz(-2.5261648) q[2];
sx q[2];
rz(-2.0912898) q[2];
sx q[2];
rz(0.25582886) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.21877737) q[1];
sx q[1];
rz(-0.73000693) q[1];
sx q[1];
rz(-3.0002712) q[1];
rz(-pi) q[2];
rz(-2.4086558) q[3];
sx q[3];
rz(-1.6896491) q[3];
sx q[3];
rz(-2.262923) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.49742302) q[2];
sx q[2];
rz(-1.3966509) q[2];
sx q[2];
rz(2.4141342) q[2];
rz(-0.71508956) q[3];
sx q[3];
rz(-0.46052027) q[3];
sx q[3];
rz(0.49811825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46399507) q[0];
sx q[0];
rz(-1.3901187) q[0];
sx q[0];
rz(-0.736262) q[0];
rz(2.5212506) q[1];
sx q[1];
rz(-0.54519975) q[1];
sx q[1];
rz(-1.5249407) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0686091) q[0];
sx q[0];
rz(-1.6260379) q[0];
sx q[0];
rz(1.5609571) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.97648804) q[2];
sx q[2];
rz(-0.39820489) q[2];
sx q[2];
rz(-1.2334241) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.273539) q[1];
sx q[1];
rz(-2.2333849) q[1];
sx q[1];
rz(2.6602547) q[1];
rz(-pi) q[2];
rz(-2.3869963) q[3];
sx q[3];
rz(-0.9661583) q[3];
sx q[3];
rz(2.5854994) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.6058495) q[2];
sx q[2];
rz(-2.3652786) q[2];
sx q[2];
rz(2.1898451) q[2];
rz(-0.4176628) q[3];
sx q[3];
rz(-2.8916736) q[3];
sx q[3];
rz(-1.1550268) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
x q[3];
rz(-pi/2) q[3];
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
rz(-0.34053892) q[0];
sx q[0];
rz(-1.8908353) q[0];
sx q[0];
rz(0.75403768) q[0];
rz(2.2484089) q[1];
sx q[1];
rz(-1.040753) q[1];
sx q[1];
rz(-2.975614) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6245859) q[0];
sx q[0];
rz(-2.4501194) q[0];
sx q[0];
rz(-0.16384478) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.31678172) q[2];
sx q[2];
rz(-2.6133279) q[2];
sx q[2];
rz(-1.007391) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.41478911) q[1];
sx q[1];
rz(-2.1926452) q[1];
sx q[1];
rz(1.0110823) q[1];
x q[2];
rz(-3.0021871) q[3];
sx q[3];
rz(-1.2981112) q[3];
sx q[3];
rz(2.8806339) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.37504998) q[2];
sx q[2];
rz(-2.0039717) q[2];
sx q[2];
rz(-0.005391187) q[2];
rz(-0.088717669) q[3];
sx q[3];
rz(-1.6088477) q[3];
sx q[3];
rz(0.79571342) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
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
rz(1.2583112) q[0];
sx q[0];
rz(-1.2092051) q[0];
sx q[0];
rz(-2.9364371) q[0];
rz(-0.908665) q[1];
sx q[1];
rz(-2.2712207) q[1];
sx q[1];
rz(-0.044513449) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9712898) q[0];
sx q[0];
rz(-1.4409541) q[0];
sx q[0];
rz(-1.4128039) q[0];
rz(-pi) q[1];
rz(-1.3088063) q[2];
sx q[2];
rz(-2.9963608) q[2];
sx q[2];
rz(-1.1373718) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.545657) q[1];
sx q[1];
rz(-2.5251198) q[1];
sx q[1];
rz(1.9484387) q[1];
rz(-pi) q[2];
rz(1.3839339) q[3];
sx q[3];
rz(-0.12178856) q[3];
sx q[3];
rz(2.361382) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-0.32555106) q[2];
sx q[2];
rz(-0.52246919) q[2];
sx q[2];
rz(-2.6204056) q[2];
rz(0.19736396) q[3];
sx q[3];
rz(-1.3857931) q[3];
sx q[3];
rz(-2.638773) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
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
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7407783) q[0];
sx q[0];
rz(-0.1939119) q[0];
sx q[0];
rz(-0.25522301) q[0];
rz(-2.9995645) q[1];
sx q[1];
rz(-1.4165001) q[1];
sx q[1];
rz(-0.35370383) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3981374) q[0];
sx q[0];
rz(-2.6530735) q[0];
sx q[0];
rz(2.3228881) q[0];
x q[1];
rz(-2.2204705) q[2];
sx q[2];
rz(-0.73474681) q[2];
sx q[2];
rz(-2.2237847) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.28359828) q[1];
sx q[1];
rz(-1.9093175) q[1];
sx q[1];
rz(-1.1095065) q[1];
rz(0.74905101) q[3];
sx q[3];
rz(-2.0045677) q[3];
sx q[3];
rz(1.6925136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.8346617) q[2];
sx q[2];
rz(-2.8262409) q[2];
sx q[2];
rz(1.6174512) q[2];
rz(0.8664242) q[3];
sx q[3];
rz(-1.6774991) q[3];
sx q[3];
rz(2.5518937) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.50188142) q[0];
sx q[0];
rz(-2.9852133) q[0];
sx q[0];
rz(1.7891275) q[0];
rz(1.9068708) q[1];
sx q[1];
rz(-2.9094978) q[1];
sx q[1];
rz(-2.629705) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0842132) q[0];
sx q[0];
rz(-1.5176) q[0];
sx q[0];
rz(1.4131111) q[0];
rz(-pi) q[1];
x q[1];
rz(1.4134679) q[2];
sx q[2];
rz(-0.81153389) q[2];
sx q[2];
rz(-3.1026156) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.7831969) q[1];
sx q[1];
rz(-1.277521) q[1];
sx q[1];
rz(0.44289987) q[1];
x q[2];
rz(-1.0591577) q[3];
sx q[3];
rz(-0.64877227) q[3];
sx q[3];
rz(1.8442475) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.0074244) q[2];
sx q[2];
rz(-1.3205426) q[2];
sx q[2];
rz(2.7443938) q[2];
rz(2.126179) q[3];
sx q[3];
rz(-2.1404603) q[3];
sx q[3];
rz(-2.5889682) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.39887244) q[0];
sx q[0];
rz(-1.9341368) q[0];
sx q[0];
rz(0.43750986) q[0];
rz(-0.73879009) q[1];
sx q[1];
rz(-2.3414325) q[1];
sx q[1];
rz(-1.4791666) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9112944) q[0];
sx q[0];
rz(-2.0574967) q[0];
sx q[0];
rz(1.507797) q[0];
rz(-pi) q[1];
x q[1];
rz(0.11585856) q[2];
sx q[2];
rz(-1.2136974) q[2];
sx q[2];
rz(3.0810205) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.4829172) q[1];
sx q[1];
rz(-1.5086996) q[1];
sx q[1];
rz(3.0351522) q[1];
rz(-2.3607632) q[3];
sx q[3];
rz(-1.025628) q[3];
sx q[3];
rz(3.0913251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.2207569) q[2];
sx q[2];
rz(-1.7831384) q[2];
sx q[2];
rz(2.9662568) q[2];
rz(-1.0026503) q[3];
sx q[3];
rz(-2.9084539) q[3];
sx q[3];
rz(-1.233915) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47068448) q[0];
sx q[0];
rz(-1.6841472) q[0];
sx q[0];
rz(2.0367784) q[0];
rz(-0.15057527) q[1];
sx q[1];
rz(-2.2655948) q[1];
sx q[1];
rz(1.2949952) q[1];
rz(2.9433123) q[2];
sx q[2];
rz(-1.5783327) q[2];
sx q[2];
rz(-0.21680149) q[2];
rz(-0.2507052) q[3];
sx q[3];
rz(-1.7964994) q[3];
sx q[3];
rz(-2.9356706) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
