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
rz(-1.3353424) q[0];
sx q[0];
rz(3.5080533) q[0];
sx q[0];
rz(8.9656497) q[0];
rz(2.4899809) q[1];
sx q[1];
rz(4.876457) q[1];
sx q[1];
rz(9.193037) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9137751) q[0];
sx q[0];
rz(-1.639606) q[0];
sx q[0];
rz(0.34389307) q[0];
x q[1];
rz(-0.094005748) q[2];
sx q[2];
rz(-2.5931907) q[2];
sx q[2];
rz(-1.3863871) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.6218839) q[1];
sx q[1];
rz(-2.1823898) q[1];
sx q[1];
rz(-2.0382463) q[1];
rz(3.0240293) q[3];
sx q[3];
rz(-2.8010578) q[3];
sx q[3];
rz(1.7528319) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.0509402) q[2];
sx q[2];
rz(-0.53477627) q[2];
sx q[2];
rz(-1.2789307) q[2];
rz(2.2517962) q[3];
sx q[3];
rz(-1.6190448) q[3];
sx q[3];
rz(0.33862996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5559674) q[0];
sx q[0];
rz(-2.2323759) q[0];
sx q[0];
rz(0.92581785) q[0];
rz(2.0766808) q[1];
sx q[1];
rz(-1.5771259) q[1];
sx q[1];
rz(0.085478641) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32498549) q[0];
sx q[0];
rz(-2.7517635) q[0];
sx q[0];
rz(0.0086602517) q[0];
rz(-pi) q[1];
rz(1.0930834) q[2];
sx q[2];
rz(-2.4730549) q[2];
sx q[2];
rz(0.82205176) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.6798415) q[1];
sx q[1];
rz(-1.0437168) q[1];
sx q[1];
rz(-2.8746469) q[1];
rz(-pi) q[2];
rz(0.68686266) q[3];
sx q[3];
rz(-1.1986889) q[3];
sx q[3];
rz(0.61249477) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-1.815879) q[2];
sx q[2];
rz(-1.4789944) q[2];
sx q[2];
rz(0.386664) q[2];
rz(-1.2169085) q[3];
sx q[3];
rz(-0.34438008) q[3];
sx q[3];
rz(2.9215422) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3306408) q[0];
sx q[0];
rz(-0.37136677) q[0];
sx q[0];
rz(2.6445342) q[0];
rz(1.0423543) q[1];
sx q[1];
rz(-0.94756871) q[1];
sx q[1];
rz(-0.29464468) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2089068) q[0];
sx q[0];
rz(-2.8021376) q[0];
sx q[0];
rz(-2.9333788) q[0];
rz(-pi) q[1];
rz(-0.23068409) q[2];
sx q[2];
rz(-0.65197021) q[2];
sx q[2];
rz(0.99562746) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-2.5897419) q[1];
sx q[1];
rz(-0.95847469) q[1];
sx q[1];
rz(2.603884) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.14400836) q[3];
sx q[3];
rz(-0.97471182) q[3];
sx q[3];
rz(-2.2726187) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-1.7872539) q[2];
sx q[2];
rz(-2.2806809) q[2];
sx q[2];
rz(-1.5798689) q[2];
rz(1.076237) q[3];
sx q[3];
rz(-1.302224) q[3];
sx q[3];
rz(2.2126183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8266066) q[0];
sx q[0];
rz(-1.5262693) q[0];
sx q[0];
rz(-0.22931799) q[0];
rz(1.6353105) q[1];
sx q[1];
rz(-1.6835667) q[1];
sx q[1];
rz(0.062072676) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6825166) q[0];
sx q[0];
rz(-2.354265) q[0];
sx q[0];
rz(0.83700257) q[0];
rz(-pi) q[1];
x q[1];
rz(1.9628635) q[2];
sx q[2];
rz(-1.3118366) q[2];
sx q[2];
rz(0.94253892) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5531299) q[1];
sx q[1];
rz(-2.3847918) q[1];
sx q[1];
rz(2.7314145) q[1];
x q[2];
rz(-2.7113879) q[3];
sx q[3];
rz(-0.22048002) q[3];
sx q[3];
rz(-1.9752432) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(0.092992358) q[2];
sx q[2];
rz(-0.91602641) q[2];
sx q[2];
rz(1.830706) q[2];
rz(0.038289573) q[3];
sx q[3];
rz(-1.3524651) q[3];
sx q[3];
rz(0.37461764) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49486092) q[0];
sx q[0];
rz(-0.72135389) q[0];
sx q[0];
rz(-0.89299655) q[0];
rz(1.0294186) q[1];
sx q[1];
rz(-2.4118377) q[1];
sx q[1];
rz(0.095887862) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.77377787) q[0];
sx q[0];
rz(-0.68499631) q[0];
sx q[0];
rz(2.2626586) q[0];
rz(-pi) q[1];
rz(-2.2402359) q[2];
sx q[2];
rz(-0.11321214) q[2];
sx q[2];
rz(2.8100546) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8687415) q[1];
sx q[1];
rz(-0.68528803) q[1];
sx q[1];
rz(2.5455695) q[1];
rz(-pi) q[2];
x q[2];
rz(0.60378243) q[3];
sx q[3];
rz(-0.60294916) q[3];
sx q[3];
rz(-0.83151885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.15478495) q[2];
sx q[2];
rz(-1.3900577) q[2];
sx q[2];
rz(2.7161157) q[2];
rz(0.42243877) q[3];
sx q[3];
rz(-1.1047624) q[3];
sx q[3];
rz(2.6384242) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.7285889) q[0];
sx q[0];
rz(-0.63218963) q[0];
sx q[0];
rz(0.11716209) q[0];
rz(1.4940184) q[1];
sx q[1];
rz(-0.90227503) q[1];
sx q[1];
rz(-1.0711627) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2881831) q[0];
sx q[0];
rz(-2.5700535) q[0];
sx q[0];
rz(2.6009212) q[0];
rz(-pi) q[1];
rz(0.44354673) q[2];
sx q[2];
rz(-2.5534592) q[2];
sx q[2];
rz(3.0556222) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.92381645) q[1];
sx q[1];
rz(-0.72715302) q[1];
sx q[1];
rz(1.9132861) q[1];
rz(-pi) q[2];
rz(1.734524) q[3];
sx q[3];
rz(-1.2571954) q[3];
sx q[3];
rz(-2.0255476) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.5693207) q[2];
sx q[2];
rz(-0.57491493) q[2];
sx q[2];
rz(-3.102109) q[2];
rz(-0.076400541) q[3];
sx q[3];
rz(-1.971784) q[3];
sx q[3];
rz(-0.45342818) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.70306784) q[0];
sx q[0];
rz(-2.330133) q[0];
sx q[0];
rz(2.1912498) q[0];
rz(1.9901265) q[1];
sx q[1];
rz(-1.5299503) q[1];
sx q[1];
rz(-1.8340402) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.81539736) q[0];
sx q[0];
rz(-1.895825) q[0];
sx q[0];
rz(-0.23120489) q[0];
rz(-pi) q[1];
x q[1];
rz(0.37613873) q[2];
sx q[2];
rz(-1.5030295) q[2];
sx q[2];
rz(0.68830953) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.89245172) q[1];
sx q[1];
rz(-2.0607407) q[1];
sx q[1];
rz(2.9016414) q[1];
x q[2];
rz(-0.96486196) q[3];
sx q[3];
rz(-1.4601225) q[3];
sx q[3];
rz(-1.0718653) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.82039708) q[2];
sx q[2];
rz(-2.4615007) q[2];
sx q[2];
rz(-0.55366984) q[2];
rz(0.6959483) q[3];
sx q[3];
rz(-1.6690994) q[3];
sx q[3];
rz(1.455201) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0221136) q[0];
sx q[0];
rz(-1.4520293) q[0];
sx q[0];
rz(-0.30690646) q[0];
rz(1.8652929) q[1];
sx q[1];
rz(-2.6752094) q[1];
sx q[1];
rz(1.9006405) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.81344714) q[0];
sx q[0];
rz(-1.7268344) q[0];
sx q[0];
rz(-2.0202083) q[0];
x q[1];
rz(0.59126258) q[2];
sx q[2];
rz(-0.68702543) q[2];
sx q[2];
rz(-2.9369773) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.7684264) q[1];
sx q[1];
rz(-1.2618999) q[1];
sx q[1];
rz(-1.7684446) q[1];
rz(-1.1697024) q[3];
sx q[3];
rz(-0.99906427) q[3];
sx q[3];
rz(2.9740262) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.79157311) q[2];
sx q[2];
rz(-2.3363523) q[2];
sx q[2];
rz(-1.8512858) q[2];
rz(2.4604515) q[3];
sx q[3];
rz(-2.2414424) q[3];
sx q[3];
rz(-1.6397938) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.27555585) q[0];
sx q[0];
rz(-2.10502) q[0];
sx q[0];
rz(1.1736897) q[0];
rz(0.58700079) q[1];
sx q[1];
rz(-2.360011) q[1];
sx q[1];
rz(3.0618844) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3088206) q[0];
sx q[0];
rz(-1.3448633) q[0];
sx q[0];
rz(1.1908047) q[0];
rz(-2.6775421) q[2];
sx q[2];
rz(-2.8949605) q[2];
sx q[2];
rz(-0.76619785) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(1.4707798) q[1];
sx q[1];
rz(-2.1683886) q[1];
sx q[1];
rz(-0.87009354) q[1];
rz(-pi) q[2];
x q[2];
rz(0.80612889) q[3];
sx q[3];
rz(-1.7793831) q[3];
sx q[3];
rz(-0.94193469) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.6012663) q[2];
sx q[2];
rz(-1.2286295) q[2];
sx q[2];
rz(-1.3339174) q[2];
rz(-1.5506844) q[3];
sx q[3];
rz(-2.1315137) q[3];
sx q[3];
rz(0.18868119) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.659336) q[0];
sx q[0];
rz(-1.5784669) q[0];
sx q[0];
rz(1.2905066) q[0];
rz(-3.105063) q[1];
sx q[1];
rz(-1.9772269) q[1];
sx q[1];
rz(2.0972924) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.8528906) q[0];
sx q[0];
rz(-1.9963309) q[0];
sx q[0];
rz(-2.8136926) q[0];
rz(-pi) q[1];
rz(-2.2453111) q[2];
sx q[2];
rz(-0.41322979) q[2];
sx q[2];
rz(2.9256224) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.83936939) q[1];
sx q[1];
rz(-1.4166792) q[1];
sx q[1];
rz(-2.3030512) q[1];
rz(-pi) q[2];
x q[2];
rz(0.32560168) q[3];
sx q[3];
rz(-0.34593098) q[3];
sx q[3];
rz(0.17445645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.2269939) q[2];
sx q[2];
rz(-2.1103766) q[2];
sx q[2];
rz(-1.8168137) q[2];
rz(-0.75657183) q[3];
sx q[3];
rz(-0.888266) q[3];
sx q[3];
rz(-0.85048401) q[3];
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
x q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.666438) q[0];
sx q[0];
rz(-1.7378687) q[0];
sx q[0];
rz(0.73874656) q[0];
rz(-1.7186164) q[1];
sx q[1];
rz(-1.1971133) q[1];
sx q[1];
rz(-2.8308629) q[1];
rz(1.995261) q[2];
sx q[2];
rz(-0.88256114) q[2];
sx q[2];
rz(0.12086856) q[2];
rz(-2.4800469) q[3];
sx q[3];
rz(-1.1894124) q[3];
sx q[3];
rz(-0.00581707) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
