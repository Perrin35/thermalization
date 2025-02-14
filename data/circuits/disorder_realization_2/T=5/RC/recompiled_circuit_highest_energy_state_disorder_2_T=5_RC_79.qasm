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
rz(1.6588563) q[0];
sx q[0];
rz(-0.98134494) q[0];
sx q[0];
rz(-2.0438097) q[0];
rz(-0.78805796) q[1];
sx q[1];
rz(-2.0907953) q[1];
sx q[1];
rz(-0.0069590574) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.46411447) q[0];
sx q[0];
rz(-2.4888746) q[0];
sx q[0];
rz(1.663289) q[0];
rz(0.9240146) q[2];
sx q[2];
rz(-1.7597464) q[2];
sx q[2];
rz(-1.3786157) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.4930058) q[1];
sx q[1];
rz(-1.3351591) q[1];
sx q[1];
rz(-3.1070903) q[1];
x q[2];
rz(1.6775292) q[3];
sx q[3];
rz(-1.9088703) q[3];
sx q[3];
rz(-2.6298912) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.9030582) q[2];
sx q[2];
rz(-1.7944585) q[2];
sx q[2];
rz(-1.4702338) q[2];
rz(3.0155731) q[3];
sx q[3];
rz(-2.6991548) q[3];
sx q[3];
rz(-0.68388763) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11494342) q[0];
sx q[0];
rz(-0.4489972) q[0];
sx q[0];
rz(-2.5057416) q[0];
rz(0.77330971) q[1];
sx q[1];
rz(-2.4644303) q[1];
sx q[1];
rz(1.6962475) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.46585402) q[0];
sx q[0];
rz(-0.93198085) q[0];
sx q[0];
rz(2.5771228) q[0];
rz(-pi) q[1];
rz(-2.3582621) q[2];
sx q[2];
rz(-1.223067) q[2];
sx q[2];
rz(1.348198) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.8938039) q[1];
sx q[1];
rz(-2.2671428) q[1];
sx q[1];
rz(-1.3651834) q[1];
rz(-pi) q[2];
rz(-1.7023193) q[3];
sx q[3];
rz(-1.5029385) q[3];
sx q[3];
rz(-2.0236286) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.14629743) q[2];
sx q[2];
rz(-1.7701021) q[2];
sx q[2];
rz(-0.13776097) q[2];
rz(-2.6855101) q[3];
sx q[3];
rz(-2.5465953) q[3];
sx q[3];
rz(-0.47541398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.544203) q[0];
sx q[0];
rz(-1.5837639) q[0];
sx q[0];
rz(-0.57918817) q[0];
rz(-3.0384565) q[1];
sx q[1];
rz(-1.8201647) q[1];
sx q[1];
rz(2.3613222) q[1];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.73585498) q[0];
sx q[0];
rz(-0.097022382) q[0];
sx q[0];
rz(2.9461622) q[0];
rz(-pi) q[1];
rz(-1.7306855) q[2];
sx q[2];
rz(-0.89582755) q[2];
sx q[2];
rz(-0.28197786) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.6578778) q[1];
sx q[1];
rz(-2.6838852) q[1];
sx q[1];
rz(1.0545516) q[1];
x q[2];
rz(-0.87522755) q[3];
sx q[3];
rz(-1.7415294) q[3];
sx q[3];
rz(0.010893498) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.4812193) q[2];
sx q[2];
rz(-1.4619091) q[2];
sx q[2];
rz(-0.090864651) q[2];
rz(-1.6615435) q[3];
sx q[3];
rz(-1.816498) q[3];
sx q[3];
rz(-0.81502325) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.018983) q[0];
sx q[0];
rz(-2.9538587) q[0];
sx q[0];
rz(-0.51938272) q[0];
rz(1.1795562) q[1];
sx q[1];
rz(-0.68125454) q[1];
sx q[1];
rz(-3.093241) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.40019401) q[0];
sx q[0];
rz(-2.1693008) q[0];
sx q[0];
rz(-0.6299751) q[0];
x q[1];
rz(1.085454) q[2];
sx q[2];
rz(-1.9159217) q[2];
sx q[2];
rz(1.3292154) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.0846588) q[1];
sx q[1];
rz(-2.6531583) q[1];
sx q[1];
rz(-0.059367511) q[1];
rz(-pi) q[2];
x q[2];
rz(0.79508852) q[3];
sx q[3];
rz(-1.3833191) q[3];
sx q[3];
rz(-2.335473) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.68253303) q[2];
sx q[2];
rz(-2.8352663) q[2];
sx q[2];
rz(1.4502067) q[2];
rz(1.1059443) q[3];
sx q[3];
rz(-1.148843) q[3];
sx q[3];
rz(-2.2753184) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.76097101) q[0];
sx q[0];
rz(-0.80046099) q[0];
sx q[0];
rz(-2.3642484) q[0];
rz(2.6457973) q[1];
sx q[1];
rz(-2.7340041) q[1];
sx q[1];
rz(0.19283238) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.042943311) q[0];
sx q[0];
rz(-2.1060364) q[0];
sx q[0];
rz(1.9283867) q[0];
rz(-pi) q[1];
x q[1];
rz(1.3309137) q[2];
sx q[2];
rz(-2.3238306) q[2];
sx q[2];
rz(0.32599923) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.3461756) q[1];
sx q[1];
rz(-2.8320791) q[1];
sx q[1];
rz(-1.8619821) q[1];
x q[2];
rz(-1.2019346) q[3];
sx q[3];
rz(-0.4274803) q[3];
sx q[3];
rz(-2.1457246) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.58482802) q[2];
sx q[2];
rz(-0.77797055) q[2];
sx q[2];
rz(0.83621109) q[2];
rz(-2.1861475) q[3];
sx q[3];
rz(-1.4232057) q[3];
sx q[3];
rz(0.53171617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8504976) q[0];
sx q[0];
rz(-1.2245155) q[0];
sx q[0];
rz(-3.1078872) q[0];
rz(-2.2233502) q[1];
sx q[1];
rz(-2.4372209) q[1];
sx q[1];
rz(-0.69923002) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3501773) q[0];
sx q[0];
rz(-0.79424131) q[0];
sx q[0];
rz(-2.3167531) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.27710931) q[2];
sx q[2];
rz(-0.19453262) q[2];
sx q[2];
rz(2.3152318) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.67656006) q[1];
sx q[1];
rz(-2.0102215) q[1];
sx q[1];
rz(1.3759173) q[1];
x q[2];
rz(-1.9101891) q[3];
sx q[3];
rz(-1.3073788) q[3];
sx q[3];
rz(-2.6201191) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0596727) q[2];
sx q[2];
rz(-2.4884188) q[2];
sx q[2];
rz(-0.095452249) q[2];
rz(-2.0898315) q[3];
sx q[3];
rz(-1.2442518) q[3];
sx q[3];
rz(-0.89422798) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8555701) q[0];
sx q[0];
rz(-0.48208553) q[0];
sx q[0];
rz(-1.1676189) q[0];
rz(-2.7883912) q[1];
sx q[1];
rz(-0.51274061) q[1];
sx q[1];
rz(-2.8599427) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1167867) q[0];
sx q[0];
rz(-1.483498) q[0];
sx q[0];
rz(2.0463819) q[0];
rz(0.64729624) q[2];
sx q[2];
rz(-1.4855097) q[2];
sx q[2];
rz(2.6814987) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.0593917) q[1];
sx q[1];
rz(-1.4363465) q[1];
sx q[1];
rz(1.7198635) q[1];
x q[2];
rz(0.48430932) q[3];
sx q[3];
rz(-2.0748027) q[3];
sx q[3];
rz(-2.8093317) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.99792751) q[2];
sx q[2];
rz(-1.9196271) q[2];
sx q[2];
rz(1.3227051) q[2];
rz(1.6392684) q[3];
sx q[3];
rz(-2.0570677) q[3];
sx q[3];
rz(3.0433906) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.1439576) q[0];
sx q[0];
rz(-1.2459545) q[0];
sx q[0];
rz(-0.27398807) q[0];
rz(2.5471089) q[1];
sx q[1];
rz(-1.7285873) q[1];
sx q[1];
rz(0.53057539) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0218567) q[0];
sx q[0];
rz(-1.6101735) q[0];
sx q[0];
rz(-1.5669797) q[0];
rz(-0.55270393) q[2];
sx q[2];
rz(-1.505902) q[2];
sx q[2];
rz(3.0187424) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1889607) q[1];
sx q[1];
rz(-2.7285693) q[1];
sx q[1];
rz(2.7160591) q[1];
x q[2];
rz(1.382393) q[3];
sx q[3];
rz(-0.72160463) q[3];
sx q[3];
rz(0.60886514) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-2.8990367) q[2];
sx q[2];
rz(-1.2071004) q[2];
sx q[2];
rz(-0.25699082) q[2];
rz(-1.1993923) q[3];
sx q[3];
rz(-1.2446087) q[3];
sx q[3];
rz(-1.6283584) q[3];
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
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5892107) q[0];
sx q[0];
rz(-2.1871545) q[0];
sx q[0];
rz(-0.5300262) q[0];
rz(1.6389305) q[1];
sx q[1];
rz(-1.9866147) q[1];
sx q[1];
rz(0.93856215) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8535271) q[0];
sx q[0];
rz(-0.99366436) q[0];
sx q[0];
rz(2.1398628) q[0];
rz(-2.9979894) q[2];
sx q[2];
rz(-1.68949) q[2];
sx q[2];
rz(2.375756) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.5903428) q[1];
sx q[1];
rz(-1.1444725) q[1];
sx q[1];
rz(-0.23386441) q[1];
x q[2];
rz(-1.3665365) q[3];
sx q[3];
rz(-2.5464749) q[3];
sx q[3];
rz(2.4581153) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-2.3878801) q[2];
sx q[2];
rz(-2.2215999) q[2];
sx q[2];
rz(-2.06125) q[2];
rz(-0.8148109) q[3];
sx q[3];
rz(-1.1956513) q[3];
sx q[3];
rz(1.4174392) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.998488) q[0];
sx q[0];
rz(-1.4838706) q[0];
sx q[0];
rz(1.0888354) q[0];
rz(3.0740652) q[1];
sx q[1];
rz(-1.2779526) q[1];
sx q[1];
rz(-2.3823104) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7831206) q[0];
sx q[0];
rz(-1.4519339) q[0];
sx q[0];
rz(-1.8274183) q[0];
rz(-pi) q[1];
rz(1.4467054) q[2];
sx q[2];
rz(-1.3465836) q[2];
sx q[2];
rz(-1.8607163) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(2.2226505) q[1];
sx q[1];
rz(-0.83178751) q[1];
sx q[1];
rz(1.2779425) q[1];
rz(-pi) q[2];
x q[2];
rz(2.5645026) q[3];
sx q[3];
rz(-0.40501696) q[3];
sx q[3];
rz(-1.7187723) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(0.33409432) q[2];
sx q[2];
rz(-1.0958025) q[2];
sx q[2];
rz(-0.57541543) q[2];
rz(-0.56257644) q[3];
sx q[3];
rz(-0.96499363) q[3];
sx q[3];
rz(1.7687198) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0811049) q[0];
sx q[0];
rz(-1.8555547) q[0];
sx q[0];
rz(-0.9077358) q[0];
rz(1.0870712) q[1];
sx q[1];
rz(-1.1320976) q[1];
sx q[1];
rz(-0.017398106) q[1];
rz(-0.74093735) q[2];
sx q[2];
rz(-0.20812427) q[2];
sx q[2];
rz(-0.32231449) q[2];
rz(2.5431332) q[3];
sx q[3];
rz(-1.2597407) q[3];
sx q[3];
rz(1.5516439) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
