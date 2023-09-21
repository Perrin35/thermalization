OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.66184008) q[0];
sx q[0];
rz(2.2979484) q[0];
sx q[0];
rz(9.2568682) q[0];
rz(-1.9703938) q[1];
sx q[1];
rz(-0.29532239) q[1];
sx q[1];
rz(3.0854316) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1355609) q[0];
sx q[0];
rz(-0.52838415) q[0];
sx q[0];
rz(2.0023268) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8499608) q[2];
sx q[2];
rz(-0.66510519) q[2];
sx q[2];
rz(-1.6601738) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.010477) q[1];
sx q[1];
rz(-1.8382204) q[1];
sx q[1];
rz(1.0130151) q[1];
x q[2];
rz(1.5264741) q[3];
sx q[3];
rz(-1.3111776) q[3];
sx q[3];
rz(2.0093902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.7636259) q[2];
sx q[2];
rz(-0.28181919) q[2];
sx q[2];
rz(-2.7089233) q[2];
rz(1.9487322) q[3];
sx q[3];
rz(-1.2377219) q[3];
sx q[3];
rz(-0.38309923) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6137961) q[0];
sx q[0];
rz(-2.6531117) q[0];
sx q[0];
rz(-1.312785) q[0];
rz(-0.20547543) q[1];
sx q[1];
rz(-0.97646362) q[1];
sx q[1];
rz(-1.9899433) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5702471) q[0];
sx q[0];
rz(-0.76128188) q[0];
sx q[0];
rz(2.0061357) q[0];
rz(-0.58354124) q[2];
sx q[2];
rz(-1.1569287) q[2];
sx q[2];
rz(-1.5577424) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.093106) q[1];
sx q[1];
rz(-1.9890607) q[1];
sx q[1];
rz(2.6622245) q[1];
rz(0.63668164) q[3];
sx q[3];
rz(-0.59906206) q[3];
sx q[3];
rz(0.77081313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.0097222086) q[2];
sx q[2];
rz(-1.6513609) q[2];
sx q[2];
rz(-2.9197664) q[2];
rz(-0.37718537) q[3];
sx q[3];
rz(-2.714034) q[3];
sx q[3];
rz(-2.3691573) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8310228) q[0];
sx q[0];
rz(-3.0492058) q[0];
sx q[0];
rz(3.1047399) q[0];
rz(0.82551461) q[1];
sx q[1];
rz(-1.3157536) q[1];
sx q[1];
rz(3.085014) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.53783137) q[0];
sx q[0];
rz(-1.0147525) q[0];
sx q[0];
rz(2.6105196) q[0];
x q[1];
rz(1.4372196) q[2];
sx q[2];
rz(-1.8250416) q[2];
sx q[2];
rz(-0.1453407) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(0.11675662) q[1];
sx q[1];
rz(-0.36326888) q[1];
sx q[1];
rz(-2.5519752) q[1];
x q[2];
rz(1.1832841) q[3];
sx q[3];
rz(-0.96252493) q[3];
sx q[3];
rz(-2.0351978) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.8686707) q[2];
sx q[2];
rz(-1.4977095) q[2];
sx q[2];
rz(-0.92612129) q[2];
rz(2.5849294) q[3];
sx q[3];
rz(-0.29354468) q[3];
sx q[3];
rz(1.042897) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.2264003) q[0];
sx q[0];
rz(-2.4214348) q[0];
sx q[0];
rz(2.3994989) q[0];
rz(2.0023951) q[1];
sx q[1];
rz(-0.4793872) q[1];
sx q[1];
rz(2.6779968) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6620561) q[0];
sx q[0];
rz(-0.89131309) q[0];
sx q[0];
rz(-2.7583073) q[0];
rz(-pi) q[1];
rz(1.5092588) q[2];
sx q[2];
rz(-1.9029641) q[2];
sx q[2];
rz(-2.5663944) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-2.3568748) q[1];
sx q[1];
rz(-0.80662913) q[1];
sx q[1];
rz(2.9033317) q[1];
x q[2];
rz(-0.47834088) q[3];
sx q[3];
rz(-1.0676749) q[3];
sx q[3];
rz(2.2174045) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(1.3670369) q[2];
sx q[2];
rz(-0.15731263) q[2];
sx q[2];
rz(-3.0920933) q[2];
rz(3.0130623) q[3];
sx q[3];
rz(-1.5481719) q[3];
sx q[3];
rz(-3.0226504) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.13609919) q[0];
sx q[0];
rz(-0.43731421) q[0];
sx q[0];
rz(-2.8438925) q[0];
rz(0.4822576) q[1];
sx q[1];
rz(-0.75459701) q[1];
sx q[1];
rz(-2.1972426) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0194861) q[0];
sx q[0];
rz(-1.6158316) q[0];
sx q[0];
rz(-1.7800063) q[0];
rz(-pi) q[1];
rz(0.40446754) q[2];
sx q[2];
rz(-2.3239115) q[2];
sx q[2];
rz(-0.58194619) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.22630616) q[1];
sx q[1];
rz(-1.5729135) q[1];
sx q[1];
rz(-1.5096942) q[1];
rz(-pi) q[2];
x q[2];
rz(0.51575757) q[3];
sx q[3];
rz(-2.6532647) q[3];
sx q[3];
rz(2.8749136) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.8828316) q[2];
sx q[2];
rz(-2.0488887) q[2];
sx q[2];
rz(0.10822254) q[2];
rz(0.0023068874) q[3];
sx q[3];
rz(-1.6121515) q[3];
sx q[3];
rz(0.32430696) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7047983) q[0];
sx q[0];
rz(-2.695485) q[0];
sx q[0];
rz(-0.58445245) q[0];
rz(2.2553518) q[1];
sx q[1];
rz(-0.61683547) q[1];
sx q[1];
rz(-3.086673) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.9960105) q[0];
sx q[0];
rz(-1.8390363) q[0];
sx q[0];
rz(3.0560188) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.6266277) q[2];
sx q[2];
rz(-2.8179114) q[2];
sx q[2];
rz(1.2602381) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.0262895) q[1];
sx q[1];
rz(-2.4541828) q[1];
sx q[1];
rz(0.61647146) q[1];
rz(-pi) q[2];
rz(-0.46338007) q[3];
sx q[3];
rz(-1.0372835) q[3];
sx q[3];
rz(-2.2802071) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.3871258) q[2];
sx q[2];
rz(-3.0209164) q[2];
sx q[2];
rz(1.0167271) q[2];
rz(2.5975442) q[3];
sx q[3];
rz(-2.7816911) q[3];
sx q[3];
rz(-1.3325161) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.465437) q[0];
sx q[0];
rz(-2.1570719) q[0];
sx q[0];
rz(2.8570535) q[0];
rz(-0.94447213) q[1];
sx q[1];
rz(-1.1962793) q[1];
sx q[1];
rz(0.91032666) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1032573) q[0];
sx q[0];
rz(-0.084780134) q[0];
sx q[0];
rz(0.8707365) q[0];
rz(-pi) q[1];
x q[1];
rz(0.34142999) q[2];
sx q[2];
rz(-1.9651946) q[2];
sx q[2];
rz(0.1614801) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0525166) q[1];
sx q[1];
rz(-2.0564582) q[1];
sx q[1];
rz(-1.7757225) q[1];
rz(-pi) q[2];
rz(2.1479285) q[3];
sx q[3];
rz(-2.2054407) q[3];
sx q[3];
rz(1.5821379) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.7515144) q[2];
sx q[2];
rz(-3.0780767) q[2];
sx q[2];
rz(2.2195623) q[2];
rz(0.56728029) q[3];
sx q[3];
rz(-1.4138979) q[3];
sx q[3];
rz(-1.0197619) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.89408016) q[0];
sx q[0];
rz(-0.67665726) q[0];
sx q[0];
rz(0.12938736) q[0];
rz(0.63240504) q[1];
sx q[1];
rz(-1.0267195) q[1];
sx q[1];
rz(-0.30050373) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43113118) q[0];
sx q[0];
rz(-0.80695242) q[0];
sx q[0];
rz(1.0459082) q[0];
x q[1];
rz(-0.95327611) q[2];
sx q[2];
rz(-2.2308908) q[2];
sx q[2];
rz(-2.2613139) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.064425163) q[1];
sx q[1];
rz(-0.38591138) q[1];
sx q[1];
rz(2.6618631) q[1];
rz(-pi) q[2];
rz(-1.237805) q[3];
sx q[3];
rz(-2.3559642) q[3];
sx q[3];
rz(-1.1803407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.5552716) q[2];
sx q[2];
rz(-0.99493146) q[2];
sx q[2];
rz(0.78197455) q[2];
rz(2.590495) q[3];
sx q[3];
rz(-1.7588153) q[3];
sx q[3];
rz(0.15792318) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.5683811) q[0];
sx q[0];
rz(-1.1567572) q[0];
sx q[0];
rz(3.0138299) q[0];
rz(-0.54221517) q[1];
sx q[1];
rz(-2.1844889) q[1];
sx q[1];
rz(2.382747) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0046783) q[0];
sx q[0];
rz(-1.5877962) q[0];
sx q[0];
rz(-1.1466115) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7636289) q[2];
sx q[2];
rz(-0.97059965) q[2];
sx q[2];
rz(2.6892975) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-0.85941852) q[1];
sx q[1];
rz(-0.71600435) q[1];
sx q[1];
rz(2.9031258) q[1];
rz(3.127029) q[3];
sx q[3];
rz(-0.90523883) q[3];
sx q[3];
rz(-1.2138106) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.0163429) q[2];
sx q[2];
rz(-1.7813851) q[2];
sx q[2];
rz(0.49003595) q[2];
rz(-1.4222493) q[3];
sx q[3];
rz(-1.9336721) q[3];
sx q[3];
rz(-1.0629883) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.35995099) q[0];
sx q[0];
rz(-2.521823) q[0];
sx q[0];
rz(-3.066257) q[0];
rz(-0.8967337) q[1];
sx q[1];
rz(-1.9672085) q[1];
sx q[1];
rz(-0.60992253) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9655351) q[0];
sx q[0];
rz(-0.81747222) q[0];
sx q[0];
rz(2.7479991) q[0];
rz(-pi) q[1];
rz(2.5209849) q[2];
sx q[2];
rz(-2.6891516) q[2];
sx q[2];
rz(-1.1021745) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.48253179) q[1];
sx q[1];
rz(-2.3606803) q[1];
sx q[1];
rz(-2.798583) q[1];
rz(-pi) q[2];
x q[2];
rz(2.0718594) q[3];
sx q[3];
rz(-0.89322972) q[3];
sx q[3];
rz(-0.60615218) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.909409) q[2];
sx q[2];
rz(-0.82513088) q[2];
sx q[2];
rz(2.4278736) q[2];
rz(-0.37832007) q[3];
sx q[3];
rz(-2.648073) q[3];
sx q[3];
rz(2.2617214) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.338035) q[0];
sx q[0];
rz(-1.9914347) q[0];
sx q[0];
rz(1.5557355) q[0];
rz(-2.4907885) q[1];
sx q[1];
rz(-1.6497859) q[1];
sx q[1];
rz(-0.12129687) q[1];
rz(1.4177633) q[2];
sx q[2];
rz(-2.9433708) q[2];
sx q[2];
rz(2.3797258) q[2];
rz(-0.25272947) q[3];
sx q[3];
rz(-2.0287632) q[3];
sx q[3];
rz(-0.91026929) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
