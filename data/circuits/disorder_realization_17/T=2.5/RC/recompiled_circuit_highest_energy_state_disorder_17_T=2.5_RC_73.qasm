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
rz(-2.1276346) q[0];
sx q[0];
rz(-1.4528217) q[0];
sx q[0];
rz(-2.5435574) q[0];
rz(1.1818089) q[1];
sx q[1];
rz(2.476517) q[1];
sx q[1];
rz(9.2718931) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6939815) q[0];
sx q[0];
rz(-1.4911528) q[0];
sx q[0];
rz(-1.682974) q[0];
rz(-pi) q[1];
x q[1];
rz(0.48798765) q[2];
sx q[2];
rz(-1.6225015) q[2];
sx q[2];
rz(-2.0590797) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.3351288) q[1];
sx q[1];
rz(-1.6963619) q[1];
sx q[1];
rz(1.7837202) q[1];
x q[2];
rz(-0.064611994) q[3];
sx q[3];
rz(-0.046053208) q[3];
sx q[3];
rz(2.8833431) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-0.87024706) q[2];
sx q[2];
rz(-0.18834867) q[2];
sx q[2];
rz(-1.9625473) q[2];
rz(-2.7315268) q[3];
sx q[3];
rz(-2.0867917) q[3];
sx q[3];
rz(-0.89455354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
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
rz(1.1882741) q[0];
sx q[0];
rz(-2.4725547) q[0];
sx q[0];
rz(0.27729312) q[0];
rz(2.9412728) q[1];
sx q[1];
rz(-0.89626139) q[1];
sx q[1];
rz(3.0576113) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0344211) q[0];
sx q[0];
rz(-2.0382529) q[0];
sx q[0];
rz(0.22308992) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7189156) q[2];
sx q[2];
rz(-2.2930055) q[2];
sx q[2];
rz(0.87004694) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(1.9878432) q[1];
sx q[1];
rz(-2.0122408) q[1];
sx q[1];
rz(2.6325339) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.6076654) q[3];
sx q[3];
rz(-1.797343) q[3];
sx q[3];
rz(-1.2372897) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(1.1818485) q[2];
sx q[2];
rz(-0.66755787) q[2];
sx q[2];
rz(2.3954771) q[2];
rz(0.6212081) q[3];
sx q[3];
rz(-0.64380232) q[3];
sx q[3];
rz(1.3889036) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(-0.32034945) q[0];
sx q[0];
rz(-0.72750434) q[0];
sx q[0];
rz(1.9387091) q[0];
rz(-2.7299643) q[1];
sx q[1];
rz(-1.2806712) q[1];
sx q[1];
rz(1.113755) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4899556) q[0];
sx q[0];
rz(-0.86511602) q[0];
sx q[0];
rz(2.5791427) q[0];
rz(-pi) q[1];
rz(2.1861548) q[2];
sx q[2];
rz(-1.9410543) q[2];
sx q[2];
rz(0.40981612) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.1542821) q[1];
sx q[1];
rz(-0.68877367) q[1];
sx q[1];
rz(1.7810526) q[1];
x q[2];
rz(2.2334072) q[3];
sx q[3];
rz(-2.326283) q[3];
sx q[3];
rz(0.96409982) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-1.0959629) q[2];
sx q[2];
rz(-2.3489504) q[2];
sx q[2];
rz(1.9795798) q[2];
rz(2.3368733) q[3];
sx q[3];
rz(-1.2730803) q[3];
sx q[3];
rz(-1.2812251) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.74748892) q[0];
sx q[0];
rz(-1.2825613) q[0];
sx q[0];
rz(1.5149186) q[0];
rz(2.6331242) q[1];
sx q[1];
rz(-2.4912806) q[1];
sx q[1];
rz(1.633684) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.29854688) q[0];
sx q[0];
rz(-0.15104476) q[0];
sx q[0];
rz(-1.5848978) q[0];
rz(-pi) q[1];
rz(-0.99782439) q[2];
sx q[2];
rz(-2.8854239) q[2];
sx q[2];
rz(-1.0333956) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(2.705756) q[1];
sx q[1];
rz(-1.8754706) q[1];
sx q[1];
rz(-0.11805822) q[1];
x q[2];
rz(-0.42801492) q[3];
sx q[3];
rz(-0.95967442) q[3];
sx q[3];
rz(-2.1111272) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-0.41370979) q[2];
sx q[2];
rz(-1.9472313) q[2];
sx q[2];
rz(-2.3940562) q[2];
rz(-1.2306635) q[3];
sx q[3];
rz(-0.80941284) q[3];
sx q[3];
rz(0.49720732) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.43664765) q[0];
sx q[0];
rz(-1.1282938) q[0];
sx q[0];
rz(1.6823912) q[0];
rz(-2.7875426) q[1];
sx q[1];
rz(-2.470128) q[1];
sx q[1];
rz(2.5669602) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1099691) q[0];
sx q[0];
rz(-1.6480443) q[0];
sx q[0];
rz(-1.5573274) q[0];
rz(-2.861205) q[2];
sx q[2];
rz(-0.78885733) q[2];
sx q[2];
rz(2.7352509) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-0.94994846) q[1];
sx q[1];
rz(-2.9086118) q[1];
sx q[1];
rz(0.22535546) q[1];
rz(-pi) q[2];
rz(-2.0784723) q[3];
sx q[3];
rz(-0.23323828) q[3];
sx q[3];
rz(2.4539029) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-2.6358801) q[2];
sx q[2];
rz(-1.9680223) q[2];
sx q[2];
rz(0.20225784) q[2];
rz(1.0328736) q[3];
sx q[3];
rz(-2.5346916) q[3];
sx q[3];
rz(0.066298299) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.1101087) q[0];
sx q[0];
rz(-2.7543572) q[0];
sx q[0];
rz(2.7701344) q[0];
rz(-0.29608852) q[1];
sx q[1];
rz(-1.5541872) q[1];
sx q[1];
rz(-2.9782226) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5510493) q[0];
sx q[0];
rz(-1.7544909) q[0];
sx q[0];
rz(0.50409533) q[0];
x q[1];
rz(-0.67739112) q[2];
sx q[2];
rz(-0.66602) q[2];
sx q[2];
rz(-0.99413727) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5886053) q[1];
sx q[1];
rz(-1.3901911) q[1];
sx q[1];
rz(0.44579472) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.3903687) q[3];
sx q[3];
rz(-2.437915) q[3];
sx q[3];
rz(2.0247519) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-2.0614193) q[2];
sx q[2];
rz(-2.8714608) q[2];
sx q[2];
rz(2.488625) q[2];
rz(0.4717007) q[3];
sx q[3];
rz(-2.3931849) q[3];
sx q[3];
rz(-0.72699839) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.93040526) q[0];
sx q[0];
rz(-2.0624332) q[0];
sx q[0];
rz(1.0373254) q[0];
rz(0.51517454) q[1];
sx q[1];
rz(-1.2778792) q[1];
sx q[1];
rz(1.3264664) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4907287) q[0];
sx q[0];
rz(-2.0959637) q[0];
sx q[0];
rz(-2.690171) q[0];
rz(0.10978077) q[2];
sx q[2];
rz(-2.7876543) q[2];
sx q[2];
rz(-2.372641) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(2.4660037) q[1];
sx q[1];
rz(-1.4741305) q[1];
sx q[1];
rz(0.76670353) q[1];
rz(-pi) q[2];
x q[2];
rz(1.2650199) q[3];
sx q[3];
rz(-1.631284) q[3];
sx q[3];
rz(-0.6636338) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-1.6933763) q[2];
sx q[2];
rz(-0.46458149) q[2];
sx q[2];
rz(0.32619897) q[2];
rz(0.97801963) q[3];
sx q[3];
rz(-1.8947442) q[3];
sx q[3];
rz(-0.090156468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
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
rz(1.6982894) q[0];
sx q[0];
rz(-2.8841618) q[0];
sx q[0];
rz(-2.8522016) q[0];
rz(3.1293213) q[1];
sx q[1];
rz(-2.8616276) q[1];
sx q[1];
rz(-1.6562921) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8809562) q[0];
sx q[0];
rz(-2.0800965) q[0];
sx q[0];
rz(-2.3348722) q[0];
x q[1];
rz(-1.5167164) q[2];
sx q[2];
rz(-1.7940494) q[2];
sx q[2];
rz(-2.8573481) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(2.0721133) q[1];
sx q[1];
rz(-0.24976191) q[1];
sx q[1];
rz(-0.62432557) q[1];
x q[2];
rz(1.1024804) q[3];
sx q[3];
rz(-1.2833929) q[3];
sx q[3];
rz(-1.772955) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.54899186) q[2];
sx q[2];
rz(-1.5697378) q[2];
sx q[2];
rz(-0.17295095) q[2];
rz(1.3744099) q[3];
sx q[3];
rz(-1.7632615) q[3];
sx q[3];
rz(1.6636728) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3953111) q[0];
sx q[0];
rz(-2.6938541) q[0];
sx q[0];
rz(2.3315499) q[0];
rz(2.7250302) q[1];
sx q[1];
rz(-0.77605334) q[1];
sx q[1];
rz(2.7966444) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.024558) q[0];
sx q[0];
rz(-0.86625615) q[0];
sx q[0];
rz(1.4092833) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.8497963) q[2];
sx q[2];
rz(-1.4115745) q[2];
sx q[2];
rz(-0.70406841) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(1.7564063) q[1];
sx q[1];
rz(-2.029772) q[1];
sx q[1];
rz(-2.7800757) q[1];
x q[2];
rz(0.79699253) q[3];
sx q[3];
rz(-2.9755962) q[3];
sx q[3];
rz(-2.5696563) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-2.714146) q[2];
sx q[2];
rz(-2.2553208) q[2];
sx q[2];
rz(0.084058849) q[2];
rz(-0.90703026) q[3];
sx q[3];
rz(-1.4182988) q[3];
sx q[3];
rz(-0.9575873) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3323988) q[0];
sx q[0];
rz(-3.0539303) q[0];
sx q[0];
rz(2.8293389) q[0];
rz(2.9088083) q[1];
sx q[1];
rz(-2.7504031) q[1];
sx q[1];
rz(-1.2871453) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2256266) q[0];
sx q[0];
rz(-0.73333987) q[0];
sx q[0];
rz(2.649113) q[0];
rz(-pi) q[1];
rz(-1.617681) q[2];
sx q[2];
rz(-2.2094036) q[2];
sx q[2];
rz(-2.435911) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.171265) q[1];
sx q[1];
rz(-2.0205775) q[1];
sx q[1];
rz(2.8697017) q[1];
rz(0.49317499) q[3];
sx q[3];
rz(-2.4552267) q[3];
sx q[3];
rz(0.92991352) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-1.1708258) q[2];
sx q[2];
rz(-1.8584741) q[2];
sx q[2];
rz(-1.4975366) q[2];
rz(-1.0981285) q[3];
sx q[3];
rz(-0.69881717) q[3];
sx q[3];
rz(-0.48925492) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.1371786) q[0];
sx q[0];
rz(-1.7902086) q[0];
sx q[0];
rz(-0.93850346) q[0];
rz(-1.3068403) q[1];
sx q[1];
rz(-0.88941457) q[1];
sx q[1];
rz(-0.79759146) q[1];
rz(2.2597792) q[2];
sx q[2];
rz(-1.4588933) q[2];
sx q[2];
rz(2.3068538) q[2];
rz(-0.14214235) q[3];
sx q[3];
rz(-2.2323043) q[3];
sx q[3];
rz(-1.4122152) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
