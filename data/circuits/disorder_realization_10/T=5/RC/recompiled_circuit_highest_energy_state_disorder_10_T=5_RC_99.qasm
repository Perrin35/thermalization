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
rz(-1.7614814) q[0];
sx q[0];
rz(-1.5505646) q[0];
sx q[0];
rz(-2.759759) q[0];
rz(2.7497357) q[1];
sx q[1];
rz(-1.5500103) q[1];
sx q[1];
rz(2.2171807) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43542433) q[0];
sx q[0];
rz(-2.2924137) q[0];
sx q[0];
rz(0.8325141) q[0];
rz(1.5906232) q[2];
sx q[2];
rz(-1.9409436) q[2];
sx q[2];
rz(-1.7006947) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.0142483) q[1];
sx q[1];
rz(-2.6963628) q[1];
sx q[1];
rz(2.8365652) q[1];
rz(0.38710164) q[3];
sx q[3];
rz(-2.0678068) q[3];
sx q[3];
rz(2.8429075) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.9367289) q[2];
sx q[2];
rz(-1.8483346) q[2];
sx q[2];
rz(1.5667685) q[2];
rz(1.1664248) q[3];
sx q[3];
rz(-2.0765442) q[3];
sx q[3];
rz(2.4488357) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4676056) q[0];
sx q[0];
rz(-2.3335712) q[0];
sx q[0];
rz(1.1760733) q[0];
rz(-3.0740956) q[1];
sx q[1];
rz(-2.564023) q[1];
sx q[1];
rz(2.8228021) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.27869895) q[0];
sx q[0];
rz(-0.79813939) q[0];
sx q[0];
rz(0.44095914) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8068878) q[2];
sx q[2];
rz(-1.1523048) q[2];
sx q[2];
rz(2.921517) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.6401854) q[1];
sx q[1];
rz(-2.1245271) q[1];
sx q[1];
rz(2.6928765) q[1];
rz(1.1652496) q[3];
sx q[3];
rz(-2.8315548) q[3];
sx q[3];
rz(-1.3998717) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.0626283) q[2];
sx q[2];
rz(-0.88125172) q[2];
sx q[2];
rz(-0.47743615) q[2];
rz(1.3580953) q[3];
sx q[3];
rz(-0.6529468) q[3];
sx q[3];
rz(-3.131598) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
rz(-0.39182144) q[0];
sx q[0];
rz(-1.3359767) q[0];
sx q[0];
rz(-1.1392449) q[0];
rz(-2.7353752) q[1];
sx q[1];
rz(-1.8832877) q[1];
sx q[1];
rz(2.4635945) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11230532) q[0];
sx q[0];
rz(-2.5480518) q[0];
sx q[0];
rz(1.4112006) q[0];
rz(-pi) q[1];
rz(1.0933206) q[2];
sx q[2];
rz(-0.68308631) q[2];
sx q[2];
rz(-0.25743279) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.19723141) q[1];
sx q[1];
rz(-1.1121096) q[1];
sx q[1];
rz(0.58802559) q[1];
x q[2];
rz(2.8077927) q[3];
sx q[3];
rz(-2.0562647) q[3];
sx q[3];
rz(0.36241058) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.71538007) q[2];
sx q[2];
rz(-1.6343626) q[2];
sx q[2];
rz(-1.0769843) q[2];
rz(-1.2601323) q[3];
sx q[3];
rz(-2.1144805) q[3];
sx q[3];
rz(2.8857359) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.28427163) q[0];
sx q[0];
rz(-1.1742598) q[0];
sx q[0];
rz(0.76048365) q[0];
rz(2.7898232) q[1];
sx q[1];
rz(-2.4302509) q[1];
sx q[1];
rz(2.9913091) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1635123) q[0];
sx q[0];
rz(-0.0038359782) q[0];
sx q[0];
rz(-2.1001966) q[0];
rz(-pi) q[1];
rz(1.1885095) q[2];
sx q[2];
rz(-1.7312382) q[2];
sx q[2];
rz(-0.40942243) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(2.6071668) q[1];
sx q[1];
rz(-2.1517506) q[1];
sx q[1];
rz(-0.39926417) q[1];
rz(-2.683879) q[3];
sx q[3];
rz(-1.8031075) q[3];
sx q[3];
rz(1.4087698) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-1.1111697) q[2];
sx q[2];
rz(-2.9141278) q[2];
sx q[2];
rz(2.5566768) q[2];
rz(-3.0758744) q[3];
sx q[3];
rz(-1.1394371) q[3];
sx q[3];
rz(-2.0871302) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86972648) q[0];
sx q[0];
rz(-1.9166742) q[0];
sx q[0];
rz(-2.4140893) q[0];
rz(1.2959405) q[1];
sx q[1];
rz(-1.0546874) q[1];
sx q[1];
rz(-2.4268699) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6520157) q[0];
sx q[0];
rz(-1.4997109) q[0];
sx q[0];
rz(2.9201304) q[0];
rz(0.3008414) q[2];
sx q[2];
rz(-1.46035) q[2];
sx q[2];
rz(3.0120603) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.7426364) q[1];
sx q[1];
rz(-2.7161437) q[1];
sx q[1];
rz(-2.0320973) q[1];
rz(-2.1470137) q[3];
sx q[3];
rz(-1.9061229) q[3];
sx q[3];
rz(-1.9063584) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.1456445) q[2];
sx q[2];
rz(-1.6137292) q[2];
sx q[2];
rz(1.145251) q[2];
rz(-0.20259419) q[3];
sx q[3];
rz(-3.0718206) q[3];
sx q[3];
rz(2.8128459) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2957434) q[0];
sx q[0];
rz(-0.29243094) q[0];
sx q[0];
rz(-0.16786815) q[0];
rz(0.56495086) q[1];
sx q[1];
rz(-1.047537) q[1];
sx q[1];
rz(1.8126743) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.560379) q[0];
sx q[0];
rz(-2.3920076) q[0];
sx q[0];
rz(-0.7313653) q[0];
x q[1];
rz(-3.0014056) q[2];
sx q[2];
rz(-1.5221688) q[2];
sx q[2];
rz(-1.0076154) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.5766474) q[1];
sx q[1];
rz(-0.35035808) q[1];
sx q[1];
rz(0.74585657) q[1];
rz(-1.8096446) q[3];
sx q[3];
rz(-1.788056) q[3];
sx q[3];
rz(-2.1715733) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.38587511) q[2];
sx q[2];
rz(-1.2473325) q[2];
sx q[2];
rz(-2.7937549) q[2];
rz(2.8268585) q[3];
sx q[3];
rz(-1.0957402) q[3];
sx q[3];
rz(-2.0033964) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.9921853) q[0];
sx q[0];
rz(-1.5062165) q[0];
sx q[0];
rz(0.95112479) q[0];
rz(-2.8490207) q[1];
sx q[1];
rz(-0.89076275) q[1];
sx q[1];
rz(2.1975885) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88144775) q[0];
sx q[0];
rz(-1.8788188) q[0];
sx q[0];
rz(0.016168895) q[0];
rz(-pi) q[1];
rz(-1.0877785) q[2];
sx q[2];
rz(-1.1406787) q[2];
sx q[2];
rz(1.5120235) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-2.2846825) q[1];
sx q[1];
rz(-1.3976685) q[1];
sx q[1];
rz(0.022927479) q[1];
x q[2];
rz(2.1066426) q[3];
sx q[3];
rz(-0.45848819) q[3];
sx q[3];
rz(1.5036086) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.49147478) q[2];
sx q[2];
rz(-2.398141) q[2];
sx q[2];
rz(-1.7000807) q[2];
rz(-0.024638351) q[3];
sx q[3];
rz(-1.325343) q[3];
sx q[3];
rz(-2.1520481) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.3330419) q[0];
sx q[0];
rz(-1.7311743) q[0];
sx q[0];
rz(-3.0035875) q[0];
rz(-2.0637312) q[1];
sx q[1];
rz(-1.4638008) q[1];
sx q[1];
rz(-0.93625751) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.9247583) q[0];
sx q[0];
rz(-2.6568597) q[0];
sx q[0];
rz(2.4527571) q[0];
rz(-0.54674863) q[2];
sx q[2];
rz(-0.90675747) q[2];
sx q[2];
rz(-0.67686096) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(3.1106134) q[1];
sx q[1];
rz(-0.18337164) q[1];
sx q[1];
rz(-1.490209) q[1];
rz(0.68558399) q[3];
sx q[3];
rz(-1.778025) q[3];
sx q[3];
rz(-0.7570467) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.1428947) q[2];
sx q[2];
rz(-1.8691209) q[2];
sx q[2];
rz(-1.8227089) q[2];
rz(-0.12254347) q[3];
sx q[3];
rz(-0.69869852) q[3];
sx q[3];
rz(-0.43012777) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.97813022) q[0];
sx q[0];
rz(-0.4438816) q[0];
sx q[0];
rz(2.7769856) q[0];
rz(1.3847903) q[1];
sx q[1];
rz(-2.2099647) q[1];
sx q[1];
rz(0.91420954) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96959463) q[0];
sx q[0];
rz(-1.9978317) q[0];
sx q[0];
rz(-2.9946505) q[0];
rz(-pi) q[1];
rz(-0.73458521) q[2];
sx q[2];
rz(-1.5365639) q[2];
sx q[2];
rz(-1.5058668) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(2.2847212) q[1];
sx q[1];
rz(-0.65509138) q[1];
sx q[1];
rz(2.9148209) q[1];
rz(0.24550415) q[3];
sx q[3];
rz(-0.55439083) q[3];
sx q[3];
rz(1.7911151) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.30283516) q[2];
sx q[2];
rz(-2.5573533) q[2];
sx q[2];
rz(-1.6287104) q[2];
rz(0.58879876) q[3];
sx q[3];
rz(-1.2062807) q[3];
sx q[3];
rz(-0.65620667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.47422472) q[0];
sx q[0];
rz(-0.56741699) q[0];
sx q[0];
rz(0.2844511) q[0];
rz(2.4868763) q[1];
sx q[1];
rz(-1.3755211) q[1];
sx q[1];
rz(0.42025748) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3271844) q[0];
sx q[0];
rz(-1.1564009) q[0];
sx q[0];
rz(0.23420686) q[0];
rz(2.3510397) q[2];
sx q[2];
rz(-1.9091064) q[2];
sx q[2];
rz(-1.6731135) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-3.0033231) q[1];
sx q[1];
rz(-1.6594634) q[1];
sx q[1];
rz(-1.6900464) q[1];
rz(-1.2154142) q[3];
sx q[3];
rz(-2.1077664) q[3];
sx q[3];
rz(-0.87600183) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.9922716) q[2];
sx q[2];
rz(-2.6565266) q[2];
sx q[2];
rz(2.1743656) q[2];
rz(1.017717) q[3];
sx q[3];
rz(-1.8571721) q[3];
sx q[3];
rz(0.73219901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.16348542) q[0];
sx q[0];
rz(-2.0757984) q[0];
sx q[0];
rz(2.1708873) q[0];
rz(-0.12374395) q[1];
sx q[1];
rz(-0.94656222) q[1];
sx q[1];
rz(1.8306517) q[1];
rz(2.8211448) q[2];
sx q[2];
rz(-1.6978227) q[2];
sx q[2];
rz(2.6200079) q[2];
rz(-1.6061892) q[3];
sx q[3];
rz(-1.8826859) q[3];
sx q[3];
rz(-1.7425698) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
