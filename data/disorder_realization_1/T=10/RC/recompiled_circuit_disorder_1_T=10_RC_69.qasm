OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.82436615) q[0];
sx q[0];
rz(5.1685652) q[0];
sx q[0];
rz(9.4246372) q[0];
rz(-1.8074942) q[1];
sx q[1];
rz(-0.9642095) q[1];
sx q[1];
rz(-1.1934086) q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0573187) q[0];
sx q[0];
rz(-2.8017375) q[0];
sx q[0];
rz(-2.1291332) q[0];
rz(2.675406) q[2];
sx q[2];
rz(-0.59980118) q[2];
sx q[2];
rz(-2.8592062) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.5692917) q[1];
sx q[1];
rz(-2.3094059) q[1];
sx q[1];
rz(0.65170793) q[1];
rz(-pi) q[2];
rz(-3.0317806) q[3];
sx q[3];
rz(-1.7870652) q[3];
sx q[3];
rz(-0.11085489) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.45941916) q[2];
sx q[2];
rz(-3.1176304) q[2];
sx q[2];
rz(1.9127282) q[2];
rz(1.4131644) q[3];
sx q[3];
rz(-2.0404405) q[3];
sx q[3];
rz(1.6536973) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6035778) q[0];
sx q[0];
rz(-1.6390272) q[0];
sx q[0];
rz(1.0128101) q[0];
rz(-3.1139328) q[1];
sx q[1];
rz(-0.67359567) q[1];
sx q[1];
rz(-1.123463) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.22966188) q[0];
sx q[0];
rz(-0.46470416) q[0];
sx q[0];
rz(-1.4421411) q[0];
rz(-1.5078817) q[2];
sx q[2];
rz(-0.79195576) q[2];
sx q[2];
rz(0.5069678) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-1.9532721) q[1];
sx q[1];
rz(-0.76474944) q[1];
sx q[1];
rz(0.79337593) q[1];
rz(-pi) q[2];
rz(-0.68617679) q[3];
sx q[3];
rz(-2.3585329) q[3];
sx q[3];
rz(-0.99883294) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(2.3479487) q[2];
sx q[2];
rz(-1.0898033) q[2];
sx q[2];
rz(-2.222555) q[2];
rz(-0.67409003) q[3];
sx q[3];
rz(-2.489311) q[3];
sx q[3];
rz(1.6154217) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8640901) q[0];
sx q[0];
rz(-0.16177495) q[0];
sx q[0];
rz(-1.2751689) q[0];
rz(-2.4480942) q[1];
sx q[1];
rz(-1.8854515) q[1];
sx q[1];
rz(2.0085874) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4988574) q[0];
sx q[0];
rz(-0.14232902) q[0];
sx q[0];
rz(-1.5940773) q[0];
rz(-pi) q[1];
rz(1.0622382) q[2];
sx q[2];
rz(-2.2937751) q[2];
sx q[2];
rz(-1.522097) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.9085711) q[1];
sx q[1];
rz(-2.1070707) q[1];
sx q[1];
rz(-2.8051393) q[1];
rz(-pi) q[2];
x q[2];
rz(0.77550019) q[3];
sx q[3];
rz(-0.79458487) q[3];
sx q[3];
rz(2.9426108) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.8901849) q[2];
sx q[2];
rz(-2.3501985) q[2];
sx q[2];
rz(1.8481002) q[2];
rz(0.039316468) q[3];
sx q[3];
rz(-1.9226363) q[3];
sx q[3];
rz(-1.2600651) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8816198) q[0];
sx q[0];
rz(-3.0631174) q[0];
sx q[0];
rz(-1.1608634) q[0];
rz(0.89598957) q[1];
sx q[1];
rz(-1.7005824) q[1];
sx q[1];
rz(3.0060351) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.7318152) q[0];
sx q[0];
rz(-1.2588132) q[0];
sx q[0];
rz(2.413495) q[0];
rz(-pi) q[1];
rz(3.0669078) q[2];
sx q[2];
rz(-1.9539781) q[2];
sx q[2];
rz(-2.243724) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.902066) q[1];
sx q[1];
rz(-2.8956928) q[1];
sx q[1];
rz(-0.86218254) q[1];
rz(2.9837708) q[3];
sx q[3];
rz(-1.2994248) q[3];
sx q[3];
rz(1.2543169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(-2.9049412) q[2];
sx q[2];
rz(-2.1950978) q[2];
sx q[2];
rz(-0.87990749) q[2];
rz(-3.0974292) q[3];
sx q[3];
rz(-1.6396089) q[3];
sx q[3];
rz(2.8529609) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0376461) q[0];
sx q[0];
rz(-2.7665311) q[0];
sx q[0];
rz(2.1283545) q[0];
rz(-3.0918616) q[1];
sx q[1];
rz(-2.2278992) q[1];
sx q[1];
rz(2.0577046) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.5729382) q[0];
sx q[0];
rz(-1.3875811) q[0];
sx q[0];
rz(-1.3230447) q[0];
rz(-1.2976228) q[2];
sx q[2];
rz(-1.3285085) q[2];
sx q[2];
rz(-0.66852409) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.10478445) q[1];
sx q[1];
rz(-2.0953062) q[1];
sx q[1];
rz(-3.0166237) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6277222) q[3];
sx q[3];
rz(-1.5268945) q[3];
sx q[3];
rz(-1.1843137) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.23285398) q[2];
sx q[2];
rz(-2.8149657) q[2];
sx q[2];
rz(0.24442913) q[2];
rz(0.43236732) q[3];
sx q[3];
rz(-1.3997388) q[3];
sx q[3];
rz(2.6385245) q[3];
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
sx q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8571092) q[0];
sx q[0];
rz(-1.720022) q[0];
sx q[0];
rz(-3.0474512) q[0];
rz(-0.17177467) q[1];
sx q[1];
rz(-1.1356907) q[1];
sx q[1];
rz(-0.89541268) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0953513) q[0];
sx q[0];
rz(-1.6086676) q[0];
sx q[0];
rz(0.34237679) q[0];
rz(-pi) q[1];
x q[1];
rz(2.5080639) q[2];
sx q[2];
rz(-1.765536) q[2];
sx q[2];
rz(-2.5031236) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-0.10455924) q[1];
sx q[1];
rz(-1.4022786) q[1];
sx q[1];
rz(3.0969572) q[1];
rz(-1.6454837) q[3];
sx q[3];
rz(-1.5684621) q[3];
sx q[3];
rz(-2.1722349) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0078997) q[2];
sx q[2];
rz(-0.40863016) q[2];
sx q[2];
rz(-2.3383979) q[2];
rz(-1.9512272) q[3];
sx q[3];
rz(-1.2322216) q[3];
sx q[3];
rz(-0.41263321) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.0728834) q[0];
sx q[0];
rz(-0.16462737) q[0];
sx q[0];
rz(-0.51914006) q[0];
rz(2.5601162) q[1];
sx q[1];
rz(-2.0362208) q[1];
sx q[1];
rz(-1.2566459) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6709267) q[0];
sx q[0];
rz(-1.5304655) q[0];
sx q[0];
rz(1.3192024) q[0];
rz(-pi) q[1];
rz(2.9314552) q[2];
sx q[2];
rz(-1.1147095) q[2];
sx q[2];
rz(-1.8535341) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.17649594) q[1];
sx q[1];
rz(-2.8526222) q[1];
sx q[1];
rz(-0.22775905) q[1];
rz(0.7092181) q[3];
sx q[3];
rz(-1.3243444) q[3];
sx q[3];
rz(2.646288) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.1533623) q[2];
sx q[2];
rz(-1.0299269) q[2];
sx q[2];
rz(1.777565) q[2];
rz(-0.91056943) q[3];
sx q[3];
rz(-1.1547337) q[3];
sx q[3];
rz(1.5301269) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
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
rz(-1.0751188) q[0];
sx q[0];
rz(-0.56448889) q[0];
sx q[0];
rz(-2.8334154) q[0];
rz(3.0691052) q[1];
sx q[1];
rz(-2.1283573) q[1];
sx q[1];
rz(2.7546308) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9403801) q[0];
sx q[0];
rz(-1.3450755) q[0];
sx q[0];
rz(-2.1696027) q[0];
rz(-pi) q[1];
rz(-0.96078028) q[2];
sx q[2];
rz(-2.3443601) q[2];
sx q[2];
rz(2.0955992) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.3124233) q[1];
sx q[1];
rz(-2.2409391) q[1];
sx q[1];
rz(1.6664684) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.49658637) q[3];
sx q[3];
rz(-1.7644617) q[3];
sx q[3];
rz(-0.95369875) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.5179634) q[2];
sx q[2];
rz(-1.364418) q[2];
sx q[2];
rz(0.44000885) q[2];
rz(-2.4258339) q[3];
sx q[3];
rz(-1.7093168) q[3];
sx q[3];
rz(-1.0796775) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
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
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0004262) q[0];
sx q[0];
rz(-2.3957802) q[0];
sx q[0];
rz(1.0986885) q[0];
rz(-0.72775841) q[1];
sx q[1];
rz(-0.37574238) q[1];
sx q[1];
rz(0.049302014) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6599777) q[0];
sx q[0];
rz(-1.0064631) q[0];
sx q[0];
rz(2.3100287) q[0];
rz(-pi) q[1];
rz(-2.2531613) q[2];
sx q[2];
rz(-1.1738136) q[2];
sx q[2];
rz(1.9156485) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.1972563) q[1];
sx q[1];
rz(-1.1068871) q[1];
sx q[1];
rz(1.0524366) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.52600577) q[3];
sx q[3];
rz(-0.81323871) q[3];
sx q[3];
rz(-0.40532743) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(3.0631642) q[2];
sx q[2];
rz(-0.59946632) q[2];
sx q[2];
rz(2.4196529) q[2];
rz(-0.94349629) q[3];
sx q[3];
rz(-2.3908581) q[3];
sx q[3];
rz(-2.8872484) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
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
rz(2.9938875) q[0];
sx q[0];
rz(-1.9839956) q[0];
sx q[0];
rz(-1.0797427) q[0];
rz(2.0823157) q[1];
sx q[1];
rz(-2.9187027) q[1];
sx q[1];
rz(-1.4019029) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.83090529) q[0];
sx q[0];
rz(-1.3584104) q[0];
sx q[0];
rz(0.12916818) q[0];
rz(-pi) q[1];
rz(-0.14342587) q[2];
sx q[2];
rz(-2.9529245) q[2];
sx q[2];
rz(2.8667237) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-0.55587308) q[1];
sx q[1];
rz(-1.8750637) q[1];
sx q[1];
rz(-0.85822206) q[1];
rz(-pi) q[2];
x q[2];
rz(0.63125061) q[3];
sx q[3];
rz(-1.4025941) q[3];
sx q[3];
rz(-1.9025308) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.4828651) q[2];
sx q[2];
rz(-1.7752703) q[2];
sx q[2];
rz(-1.6213017) q[2];
rz(2.5907717) q[3];
sx q[3];
rz(-2.3362624) q[3];
sx q[3];
rz(-0.6974535) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14810066) q[0];
sx q[0];
rz(-1.8363331) q[0];
sx q[0];
rz(1.6114417) q[0];
rz(-2.2254754) q[1];
sx q[1];
rz(-0.59090186) q[1];
sx q[1];
rz(-0.59060243) q[1];
rz(-2.3475636) q[2];
sx q[2];
rz(-2.2326438) q[2];
sx q[2];
rz(2.2868962) q[2];
rz(1.6021452) q[3];
sx q[3];
rz(-2.6242704) q[3];
sx q[3];
rz(-0.27192413) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];