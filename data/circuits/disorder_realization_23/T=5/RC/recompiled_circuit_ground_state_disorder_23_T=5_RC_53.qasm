OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(0.57103676) q[0];
sx q[0];
rz(3.8138226) q[0];
sx q[0];
rz(11.091118) q[0];
rz(0.9530468) q[1];
sx q[1];
rz(-2.9763728) q[1];
sx q[1];
rz(1.2686977) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8333913) q[0];
sx q[0];
rz(-1.9692288) q[0];
sx q[0];
rz(-1.3223668) q[0];
rz(-0.95999865) q[2];
sx q[2];
rz(-0.75163254) q[2];
sx q[2];
rz(-0.014873504) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(0.70763146) q[1];
sx q[1];
rz(-0.58506706) q[1];
sx q[1];
rz(1.8088887) q[1];
rz(-pi) q[2];
rz(2.8976909) q[3];
sx q[3];
rz(-2.3696405) q[3];
sx q[3];
rz(0.011079196) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.1573726) q[2];
sx q[2];
rz(-0.34175384) q[2];
sx q[2];
rz(2.3483707) q[2];
rz(-0.67304099) q[3];
sx q[3];
rz(-2.0561736) q[3];
sx q[3];
rz(1.8719155) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.15758812) q[0];
sx q[0];
rz(-2.3389811) q[0];
sx q[0];
rz(2.5322835) q[0];
rz(2.9318103) q[1];
sx q[1];
rz(-2.3502626) q[1];
sx q[1];
rz(0.9598859) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.083733746) q[0];
sx q[0];
rz(-1.4918033) q[0];
sx q[0];
rz(1.6154438) q[0];
rz(1.0823971) q[2];
sx q[2];
rz(-1.5090183) q[2];
sx q[2];
rz(0.90324963) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-2.7916225) q[1];
sx q[1];
rz(-2.4020139) q[1];
sx q[1];
rz(1.6639568) q[1];
rz(-pi) q[2];
x q[2];
rz(1.0632564) q[3];
sx q[3];
rz(-1.9572658) q[3];
sx q[3];
rz(1.8563351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.7858872) q[2];
sx q[2];
rz(-2.3023119) q[2];
sx q[2];
rz(-2.6500224) q[2];
rz(3.135904) q[3];
sx q[3];
rz(-1.0892884) q[3];
sx q[3];
rz(2.6277241) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.095057644) q[0];
sx q[0];
rz(-2.6061366) q[0];
sx q[0];
rz(0.74263483) q[0];
rz(0.62478089) q[1];
sx q[1];
rz(-0.94015986) q[1];
sx q[1];
rz(2.0754441) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2258573) q[0];
sx q[0];
rz(-1.7122147) q[0];
sx q[0];
rz(-1.7959602) q[0];
rz(-pi) q[1];
x q[1];
rz(1.0451308) q[2];
sx q[2];
rz(-2.774775) q[2];
sx q[2];
rz(-0.17772533) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.448145) q[1];
sx q[1];
rz(-1.9733917) q[1];
sx q[1];
rz(3.0047396) q[1];
x q[2];
rz(0.2188628) q[3];
sx q[3];
rz(-1.0893679) q[3];
sx q[3];
rz(-2.4784513) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.2523969) q[2];
sx q[2];
rz(-2.9434581) q[2];
sx q[2];
rz(0.55257094) q[2];
rz(0.50734723) q[3];
sx q[3];
rz(-1.0891958) q[3];
sx q[3];
rz(-0.56940091) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.4813389) q[0];
sx q[0];
rz(-2.5642671) q[0];
sx q[0];
rz(-1.8950155) q[0];
rz(1.4056816) q[1];
sx q[1];
rz(-2.3697) q[1];
sx q[1];
rz(2.4235922) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.8193977) q[0];
sx q[0];
rz(-0.85032636) q[0];
sx q[0];
rz(-2.2706991) q[0];
rz(-1.0077072) q[2];
sx q[2];
rz(-1.3674842) q[2];
sx q[2];
rz(2.2846534) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-0.79149063) q[1];
sx q[1];
rz(-1.7161421) q[1];
sx q[1];
rz(-0.75993698) q[1];
rz(0.15851373) q[3];
sx q[3];
rz(-1.6974546) q[3];
sx q[3];
rz(-1.4617006) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(1.8743073) q[2];
sx q[2];
rz(-2.3882046) q[2];
sx q[2];
rz(-0.70449746) q[2];
rz(2.7590175) q[3];
sx q[3];
rz(-1.7035328) q[3];
sx q[3];
rz(-2.2883435) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.48069561) q[0];
sx q[0];
rz(-0.040204164) q[0];
sx q[0];
rz(-0.30886343) q[0];
rz(-0.94912306) q[1];
sx q[1];
rz(-2.821065) q[1];
sx q[1];
rz(2.5905051) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.26924414) q[0];
sx q[0];
rz(-1.9595032) q[0];
sx q[0];
rz(1.7354911) q[0];
rz(-1.6597719) q[2];
sx q[2];
rz(-2.4818261) q[2];
sx q[2];
rz(-0.24487409) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(-2.2998553) q[1];
sx q[1];
rz(-1.4159214) q[1];
sx q[1];
rz(2.9938404) q[1];
rz(2.3651028) q[3];
sx q[3];
rz(-1.6740419) q[3];
sx q[3];
rz(-0.50738996) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(2.4264195) q[2];
sx q[2];
rz(-0.74843633) q[2];
sx q[2];
rz(0.61881649) q[2];
rz(-0.62301451) q[3];
sx q[3];
rz(-0.84191936) q[3];
sx q[3];
rz(-2.714341) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.1022559) q[0];
sx q[0];
rz(-2.4781041) q[0];
sx q[0];
rz(-2.6797507) q[0];
rz(2.6448008) q[1];
sx q[1];
rz(-1.3235612) q[1];
sx q[1];
rz(-1.7740446) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.38765466) q[0];
sx q[0];
rz(-0.68635041) q[0];
sx q[0];
rz(-1.4453056) q[0];
rz(-2.9890247) q[2];
sx q[2];
rz(-1.0786453) q[2];
sx q[2];
rz(-2.9967665) q[2];
rz(pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.16713472) q[1];
sx q[1];
rz(-0.30213812) q[1];
sx q[1];
rz(-1.4767417) q[1];
rz(-pi) q[2];
rz(2.9799906) q[3];
sx q[3];
rz(-2.186612) q[3];
sx q[3];
rz(3.1231669) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.10609145) q[2];
sx q[2];
rz(-1.1320628) q[2];
sx q[2];
rz(0.79037017) q[2];
rz(-0.60449374) q[3];
sx q[3];
rz(-1.0249745) q[3];
sx q[3];
rz(-2.7214971) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.4496434) q[0];
sx q[0];
rz(-0.88570166) q[0];
sx q[0];
rz(0.41437909) q[0];
rz(-2.5212133) q[1];
sx q[1];
rz(-1.9199771) q[1];
sx q[1];
rz(1.5588123) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-3.1315077) q[0];
sx q[0];
rz(-1.4580618) q[0];
sx q[0];
rz(-1.6555863) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.10509872) q[2];
sx q[2];
rz(-1.5529648) q[2];
sx q[2];
rz(0.43184973) q[2];
rz(-pi) q[3];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-0.95629809) q[1];
sx q[1];
rz(-1.6285588) q[1];
sx q[1];
rz(3.1321146) q[1];
rz(0.15191468) q[3];
sx q[3];
rz(-1.9524116) q[3];
sx q[3];
rz(-1.8982062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.9862426) q[2];
sx q[2];
rz(-2.4537931) q[2];
sx q[2];
rz(-2.9637994) q[2];
rz(2.6627461) q[3];
sx q[3];
rz(-1.1136473) q[3];
sx q[3];
rz(3.0732885) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
x q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.14684045) q[0];
sx q[0];
rz(-0.9786334) q[0];
sx q[0];
rz(0.11257182) q[0];
rz(2.5162137) q[1];
sx q[1];
rz(-1.1275147) q[1];
sx q[1];
rz(2.2523527) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.205394) q[0];
sx q[0];
rz(-1.5279211) q[0];
sx q[0];
rz(-1.6498736) q[0];
rz(2.8761112) q[2];
sx q[2];
rz(-1.4944139) q[2];
sx q[2];
rz(-0.35778174) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.1804643) q[1];
sx q[1];
rz(-2.5601009) q[1];
sx q[1];
rz(1.1322458) q[1];
rz(-pi) q[2];
rz(3.0744138) q[3];
sx q[3];
rz(-1.805873) q[3];
sx q[3];
rz(-2.377411) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-1.858295) q[2];
sx q[2];
rz(-0.12792835) q[2];
sx q[2];
rz(0.27321401) q[2];
rz(2.9122399) q[3];
sx q[3];
rz(-1.5867686) q[3];
sx q[3];
rz(-1.8730414) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.6943618) q[0];
sx q[0];
rz(-0.87089592) q[0];
sx q[0];
rz(-0.91621512) q[0];
rz(-2.9592196) q[1];
sx q[1];
rz(-0.20621754) q[1];
sx q[1];
rz(1.1590385) q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(3.0908546) q[0];
sx q[0];
rz(-1.3373475) q[0];
sx q[0];
rz(-1.3755193) q[0];
rz(3.0460814) q[2];
sx q[2];
rz(-1.6217124) q[2];
sx q[2];
rz(-1.7820304) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-2.7024544) q[1];
sx q[1];
rz(-1.7299486) q[1];
sx q[1];
rz(-1.8887146) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.69144122) q[3];
sx q[3];
rz(-2.2855298) q[3];
sx q[3];
rz(-0.5428398) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.75659043) q[2];
sx q[2];
rz(-0.80005163) q[2];
sx q[2];
rz(0.315256) q[2];
rz(-2.5473525) q[3];
sx q[3];
rz(-2.2599594) q[3];
sx q[3];
rz(-2.5074904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.88157982) q[0];
sx q[0];
rz(-0.6686815) q[0];
sx q[0];
rz(0.15923937) q[0];
rz(-2.0533994) q[1];
sx q[1];
rz(-0.91251487) q[1];
sx q[1];
rz(-0.27096567) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.78346398) q[0];
sx q[0];
rz(-2.6081351) q[0];
sx q[0];
rz(2.0808043) q[0];
x q[1];
rz(1.4681533) q[2];
sx q[2];
rz(-2.2410975) q[2];
sx q[2];
rz(1.070529) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(1.7587563) q[1];
sx q[1];
rz(-2.3985632) q[1];
sx q[1];
rz(-0.50937517) q[1];
rz(2.5529717) q[3];
sx q[3];
rz(-0.56097066) q[3];
sx q[3];
rz(-0.64631337) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(1.3403885) q[2];
sx q[2];
rz(-0.77216721) q[2];
sx q[2];
rz(-2.8188952) q[2];
rz(-0.05803756) q[3];
sx q[3];
rz(-2.3400584) q[3];
sx q[3];
rz(2.8431852) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4651481) q[0];
sx q[0];
rz(-1.6315176) q[0];
sx q[0];
rz(2.3254707) q[0];
rz(-0.25794087) q[1];
sx q[1];
rz(-1.1358658) q[1];
sx q[1];
rz(1.6075016) q[1];
rz(1.8564687) q[2];
sx q[2];
rz(-1.8485879) q[2];
sx q[2];
rz(-1.585203) q[2];
rz(2.4856223) q[3];
sx q[3];
rz(-1.9269939) q[3];
sx q[3];
rz(-3.1300598) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
