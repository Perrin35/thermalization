OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.82233959) q[0];
sx q[0];
rz(5.3263721) q[0];
sx q[0];
rz(10.858067) q[0];
rz(-0.23694555) q[1];
sx q[1];
rz(4.3720923) q[1];
sx q[1];
rz(11.421539) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5231397) q[0];
sx q[0];
rz(-1.8123527) q[0];
sx q[0];
rz(3.0735344) q[0];
rz(-pi) q[1];
rz(0.28540622) q[2];
sx q[2];
rz(-1.1661582) q[2];
sx q[2];
rz(-1.0813431) q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.45935985) q[1];
sx q[1];
rz(-2.0233005) q[1];
sx q[1];
rz(0.020999055) q[1];
rz(2.5979795) q[3];
sx q[3];
rz(-1.8279148) q[3];
sx q[3];
rz(1.6888113) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-3.0248489) q[2];
sx q[2];
rz(-1.3329093) q[2];
sx q[2];
rz(1.9072745) q[2];
rz(-2.0862789) q[3];
sx q[3];
rz(-2.0827259) q[3];
sx q[3];
rz(-1.1528667) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
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
rz(2.9556483) q[0];
sx q[0];
rz(-1.1971104) q[0];
sx q[0];
rz(-0.82988513) q[0];
rz(0.25289598) q[1];
sx q[1];
rz(-0.77650944) q[1];
sx q[1];
rz(-0.84709644) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0704437) q[0];
sx q[0];
rz(-2.0111472) q[0];
sx q[0];
rz(1.1308934) q[0];
rz(-pi) q[1];
rz(2.8424938) q[2];
sx q[2];
rz(-2.3397589) q[2];
sx q[2];
rz(2.9849844) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.4699128) q[1];
sx q[1];
rz(-0.54930733) q[1];
sx q[1];
rz(3.0104907) q[1];
x q[2];
rz(-0.91244016) q[3];
sx q[3];
rz(-2.6077301) q[3];
sx q[3];
rz(0.57953366) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.117924) q[2];
sx q[2];
rz(-2.7682722) q[2];
sx q[2];
rz(-0.71933293) q[2];
rz(-1.2891278) q[3];
sx q[3];
rz(-2.9057861) q[3];
sx q[3];
rz(1.6524564) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8711202) q[0];
sx q[0];
rz(-1.8914762) q[0];
sx q[0];
rz(2.5701994) q[0];
rz(-2.1748623) q[1];
sx q[1];
rz(-1.2582018) q[1];
sx q[1];
rz(2.6142696) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.2741094) q[0];
sx q[0];
rz(-1.2615146) q[0];
sx q[0];
rz(1.563619) q[0];
x q[1];
rz(2.7408319) q[2];
sx q[2];
rz(-0.71360613) q[2];
sx q[2];
rz(-2.1378627) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-2.5273683) q[1];
sx q[1];
rz(-0.80263153) q[1];
sx q[1];
rz(-0.36348344) q[1];
rz(0.35964386) q[3];
sx q[3];
rz(-2.375964) q[3];
sx q[3];
rz(-1.5776724) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-0.30806914) q[2];
sx q[2];
rz(-1.4493194) q[2];
sx q[2];
rz(-0.72675881) q[2];
rz(-2.1440992) q[3];
sx q[3];
rz(-0.62524978) q[3];
sx q[3];
rz(1.5184901) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8961287) q[0];
sx q[0];
rz(-1.6587057) q[0];
sx q[0];
rz(-2.1380651) q[0];
rz(3.1009122) q[1];
sx q[1];
rz(-1.129312) q[1];
sx q[1];
rz(-0.8262659) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.47047939) q[0];
sx q[0];
rz(-2.2320036) q[0];
sx q[0];
rz(0.83755042) q[0];
rz(-pi) q[1];
rz(0.063886558) q[2];
sx q[2];
rz(-2.4404844) q[2];
sx q[2];
rz(-1.0039767) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.1288209) q[1];
sx q[1];
rz(-1.2342617) q[1];
sx q[1];
rz(-0.080288447) q[1];
rz(-pi) q[2];
rz(1.9598947) q[3];
sx q[3];
rz(-2.1660921) q[3];
sx q[3];
rz(0.94483313) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.56759175) q[2];
sx q[2];
rz(-1.2446128) q[2];
sx q[2];
rz(0.91439247) q[2];
rz(-3.0363723) q[3];
sx q[3];
rz(-1.6243694) q[3];
sx q[3];
rz(0.58925327) q[3];
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
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88702622) q[0];
sx q[0];
rz(-1.6396602) q[0];
sx q[0];
rz(-2.0078833) q[0];
rz(3.1359613) q[1];
sx q[1];
rz(-1.7040323) q[1];
sx q[1];
rz(2.5240135) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78248258) q[0];
sx q[0];
rz(-1.6164653) q[0];
sx q[0];
rz(-0.029687667) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.4992906) q[2];
sx q[2];
rz(-1.2080492) q[2];
sx q[2];
rz(-1.158266) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.6286271) q[1];
sx q[1];
rz(-1.291853) q[1];
sx q[1];
rz(-2.1668424) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.17554749) q[3];
sx q[3];
rz(-1.5834337) q[3];
sx q[3];
rz(0.73568425) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.3474943) q[2];
sx q[2];
rz(-1.8811052) q[2];
sx q[2];
rz(0.23362544) q[2];
rz(2.1485093) q[3];
sx q[3];
rz(-0.94998327) q[3];
sx q[3];
rz(2.5036507) q[3];
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
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.1722906) q[0];
sx q[0];
rz(-0.06047051) q[0];
sx q[0];
rz(3.0517975) q[0];
rz(-2.2604997) q[1];
sx q[1];
rz(-0.52905622) q[1];
sx q[1];
rz(0.18009137) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3072309) q[0];
sx q[0];
rz(-2.0198856) q[0];
sx q[0];
rz(2.2063072) q[0];
rz(-pi) q[1];
x q[1];
rz(1.8413576) q[2];
sx q[2];
rz(-2.5205043) q[2];
sx q[2];
rz(-1.6294711) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(3.0568911) q[1];
sx q[1];
rz(-1.4086205) q[1];
sx q[1];
rz(-0.2918891) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.6052386) q[3];
sx q[3];
rz(-2.648571) q[3];
sx q[3];
rz(0.29355129) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(2.883541) q[2];
sx q[2];
rz(-1.5641314) q[2];
sx q[2];
rz(1.9656666) q[2];
rz(-0.073143395) q[3];
sx q[3];
rz(-2.4172343) q[3];
sx q[3];
rz(-0.49889645) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.183855) q[0];
sx q[0];
rz(-0.65559214) q[0];
sx q[0];
rz(1.7234329) q[0];
rz(-1.9765967) q[1];
sx q[1];
rz(-2.5823451) q[1];
sx q[1];
rz(3.1013536) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.2049853) q[0];
sx q[0];
rz(-1.4892704) q[0];
sx q[0];
rz(-1.3315014) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.1646541) q[2];
sx q[2];
rz(-1.5283661) q[2];
sx q[2];
rz(-0.050780642) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(0.38899598) q[1];
sx q[1];
rz(-2.1850359) q[1];
sx q[1];
rz(2.9700301) q[1];
rz(0.76544806) q[3];
sx q[3];
rz(-2.7141889) q[3];
sx q[3];
rz(-3.1020853) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.13850257) q[2];
sx q[2];
rz(-0.78531992) q[2];
sx q[2];
rz(2.9677532) q[2];
rz(-1.396817) q[3];
sx q[3];
rz(-1.6766179) q[3];
sx q[3];
rz(2.2824536) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1180856) q[0];
sx q[0];
rz(-0.21723391) q[0];
sx q[0];
rz(-1.7370976) q[0];
rz(-1.9346168) q[1];
sx q[1];
rz(-1.6563621) q[1];
sx q[1];
rz(1.6361902) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.36639402) q[0];
sx q[0];
rz(-1.9885855) q[0];
sx q[0];
rz(-1.9154857) q[0];
rz(0.82377388) q[2];
sx q[2];
rz(-2.0588377) q[2];
sx q[2];
rz(1.4866231) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.1464403) q[1];
sx q[1];
rz(-2.4447828) q[1];
sx q[1];
rz(-3.0909096) q[1];
rz(2.0309615) q[3];
sx q[3];
rz(-0.65350973) q[3];
sx q[3];
rz(-1.6516407) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.3796842) q[2];
sx q[2];
rz(-0.48019335) q[2];
sx q[2];
rz(-2.9910679) q[2];
rz(-1.6020417) q[3];
sx q[3];
rz(-1.9463041) q[3];
sx q[3];
rz(-2.5185744) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
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
rz(-1.4432916) q[0];
sx q[0];
rz(-2.0069831) q[0];
sx q[0];
rz(0.64569965) q[0];
rz(2.6422016) q[1];
sx q[1];
rz(-1.4285409) q[1];
sx q[1];
rz(1.8766778) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0255233) q[0];
sx q[0];
rz(-2.2795296) q[0];
sx q[0];
rz(-0.52225964) q[0];
x q[1];
rz(-0.64847704) q[2];
sx q[2];
rz(-2.4816374) q[2];
sx q[2];
rz(-2.6476268) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.7449194) q[1];
sx q[1];
rz(-1.5349689) q[1];
sx q[1];
rz(-0.14573914) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.42322741) q[3];
sx q[3];
rz(-1.3864577) q[3];
sx q[3];
rz(1.9870027) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.50679961) q[2];
sx q[2];
rz(-0.96968499) q[2];
sx q[2];
rz(0.52465049) q[2];
rz(-2.5850885) q[3];
sx q[3];
rz(-1.6289214) q[3];
sx q[3];
rz(2.7867253) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(-1.6181347) q[0];
sx q[0];
rz(-0.83704656) q[0];
sx q[0];
rz(0.24542228) q[0];
rz(-0.55631176) q[1];
sx q[1];
rz(-1.5309155) q[1];
sx q[1];
rz(1.4642749) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3905555) q[0];
sx q[0];
rz(-1.7883736) q[0];
sx q[0];
rz(1.6839954) q[0];
rz(-0.22428959) q[2];
sx q[2];
rz(-0.99594342) q[2];
sx q[2];
rz(-1.6155417) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.7732685) q[1];
sx q[1];
rz(-1.3001633) q[1];
sx q[1];
rz(2.2776105) q[1];
rz(-pi) q[2];
x q[2];
rz(0.085300307) q[3];
sx q[3];
rz(-1.395426) q[3];
sx q[3];
rz(2.8773957) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.2422553) q[2];
sx q[2];
rz(-1.9031886) q[2];
sx q[2];
rz(2.396092) q[2];
rz(1.4238822) q[3];
sx q[3];
rz(-0.29885492) q[3];
sx q[3];
rz(2.0132813) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8650919) q[0];
sx q[0];
rz(-1.4958953) q[0];
sx q[0];
rz(-2.4687742) q[0];
rz(-1.8854234) q[1];
sx q[1];
rz(-0.80614631) q[1];
sx q[1];
rz(2.0731906) q[1];
rz(1.0829265) q[2];
sx q[2];
rz(-2.2219873) q[2];
sx q[2];
rz(1.7359003) q[2];
rz(-1.9477378) q[3];
sx q[3];
rz(-2.1115163) q[3];
sx q[3];
rz(1.9797309) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
