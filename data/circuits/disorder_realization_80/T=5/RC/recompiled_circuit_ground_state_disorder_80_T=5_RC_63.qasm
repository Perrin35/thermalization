OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.4898981) q[0];
sx q[0];
rz(-2.2597921) q[0];
sx q[0];
rz(-2.9969126) q[0];
rz(0.41809234) q[1];
sx q[1];
rz(-0.10377181) q[1];
sx q[1];
rz(-1.6296847) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.9448175) q[0];
sx q[0];
rz(-0.51048764) q[0];
sx q[0];
rz(-1.8148242) q[0];
rz(1.51654) q[2];
sx q[2];
rz(-0.83803229) q[2];
sx q[2];
rz(0.65544477) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-0.92951951) q[1];
sx q[1];
rz(-1.9572502) q[1];
sx q[1];
rz(-0.35915908) q[1];
rz(-pi) q[2];
x q[2];
rz(2.469796) q[3];
sx q[3];
rz(-0.62272391) q[3];
sx q[3];
rz(-2.0178362) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(2.8251553) q[2];
sx q[2];
rz(-1.8708159) q[2];
sx q[2];
rz(-0.6553418) q[2];
rz(1.7573028) q[3];
sx q[3];
rz(-0.96450788) q[3];
sx q[3];
rz(0.82096076) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
rz(-pi) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.4493988) q[0];
sx q[0];
rz(-1.7658424) q[0];
sx q[0];
rz(-2.6334515) q[0];
rz(1.3805768) q[1];
sx q[1];
rz(-2.5602129) q[1];
sx q[1];
rz(-1.2095721) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.094645) q[0];
sx q[0];
rz(-2.1535465) q[0];
sx q[0];
rz(0.21999448) q[0];
rz(-pi) q[1];
x q[1];
rz(3.0879966) q[2];
sx q[2];
rz(-1.5531248) q[2];
sx q[2];
rz(1.1143717) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.75825071) q[1];
sx q[1];
rz(-1.1102601) q[1];
sx q[1];
rz(1.7136445) q[1];
x q[2];
rz(2.1048412) q[3];
sx q[3];
rz(-0.67592127) q[3];
sx q[3];
rz(0.40646857) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.2035344) q[2];
sx q[2];
rz(-1.396023) q[2];
sx q[2];
rz(2.5246942) q[2];
rz(-1.6054224) q[3];
sx q[3];
rz(-2.0974396) q[3];
sx q[3];
rz(-0.49083403) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9817552) q[0];
sx q[0];
rz(-1.497739) q[0];
sx q[0];
rz(0.72506654) q[0];
rz(-1.8720576) q[1];
sx q[1];
rz(-0.67765647) q[1];
sx q[1];
rz(0.95917541) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0051668) q[0];
sx q[0];
rz(-2.040103) q[0];
sx q[0];
rz(-0.43715663) q[0];
x q[1];
rz(-2.3940635) q[2];
sx q[2];
rz(-0.5420712) q[2];
sx q[2];
rz(-0.67491787) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-0.19634464) q[1];
sx q[1];
rz(-1.7540534) q[1];
sx q[1];
rz(-3.0250579) q[1];
rz(-pi) q[2];
x q[2];
rz(2.3183075) q[3];
sx q[3];
rz(-0.24634493) q[3];
sx q[3];
rz(2.5096517) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(1.8150669) q[2];
sx q[2];
rz(-1.8748137) q[2];
sx q[2];
rz(-0.25700021) q[2];
rz(2.0604996) q[3];
sx q[3];
rz(-2.6205781) q[3];
sx q[3];
rz(1.5628373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
sx q[3];
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
rz(2.0576039) q[0];
sx q[0];
rz(-0.27844089) q[0];
sx q[0];
rz(1.9637015) q[0];
rz(-2.9914757) q[1];
sx q[1];
rz(-0.53214407) q[1];
sx q[1];
rz(1.5163126) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.7895607) q[0];
sx q[0];
rz(-1.2403249) q[0];
sx q[0];
rz(-0.29833692) q[0];
x q[1];
rz(0.37090918) q[2];
sx q[2];
rz(-2.264655) q[2];
sx q[2];
rz(0.81028509) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.5387568) q[1];
sx q[1];
rz(-1.2456166) q[1];
sx q[1];
rz(-1.4657337) q[1];
rz(-pi) q[2];
rz(-2.4978906) q[3];
sx q[3];
rz(-1.7233522) q[3];
sx q[3];
rz(1.8495454) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(0.84245044) q[2];
sx q[2];
rz(-2.5717042) q[2];
sx q[2];
rz(-2.2184856) q[2];
rz(-2.182377) q[3];
sx q[3];
rz(-2.1289289) q[3];
sx q[3];
rz(-2.925351) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.32090309) q[0];
sx q[0];
rz(-2.237759) q[0];
sx q[0];
rz(2.2586816) q[0];
rz(-1.2954905) q[1];
sx q[1];
rz(-2.3675282) q[1];
sx q[1];
rz(2.6466218) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.2705602) q[0];
sx q[0];
rz(-2.705276) q[0];
sx q[0];
rz(0.83026921) q[0];
rz(-pi) q[1];
x q[1];
rz(0.70303838) q[2];
sx q[2];
rz(-1.1482757) q[2];
sx q[2];
rz(0.51681821) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(0.40303142) q[1];
sx q[1];
rz(-2.237169) q[1];
sx q[1];
rz(1.6825466) q[1];
rz(-pi) q[2];
rz(-2.0932156) q[3];
sx q[3];
rz(-1.8847147) q[3];
sx q[3];
rz(-2.3274904) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(1.8098658) q[2];
sx q[2];
rz(-0.11955424) q[2];
sx q[2];
rz(-1.8734107) q[2];
rz(-3.0590893) q[3];
sx q[3];
rz(-1.637633) q[3];
sx q[3];
rz(-0.65703195) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
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
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6562011) q[0];
sx q[0];
rz(-2.9347561) q[0];
sx q[0];
rz(1.50151) q[0];
rz(2.0791176) q[1];
sx q[1];
rz(-1.9891918) q[1];
sx q[1];
rz(-1.9662439) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0153246) q[0];
sx q[0];
rz(-1.7909174) q[0];
sx q[0];
rz(-0.31173978) q[0];
rz(-pi) q[1];
rz(-1.2771036) q[2];
sx q[2];
rz(-2.2429464) q[2];
sx q[2];
rz(-1.9543242) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-2.8379587) q[1];
sx q[1];
rz(-2.8629747) q[1];
sx q[1];
rz(0.4391123) q[1];
rz(-pi) q[2];
rz(1.8683942) q[3];
sx q[3];
rz(-2.3087325) q[3];
sx q[3];
rz(1.1296062) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(0.33092734) q[2];
sx q[2];
rz(-1.8570447) q[2];
sx q[2];
rz(2.1539099) q[2];
rz(-2.2641613) q[3];
sx q[3];
rz(-2.8743447) q[3];
sx q[3];
rz(-2.4115244) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[3];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1838609) q[0];
sx q[0];
rz(-1.6022302) q[0];
sx q[0];
rz(3.0027332) q[0];
rz(1.9271556) q[1];
sx q[1];
rz(-1.5077533) q[1];
sx q[1];
rz(0.81659281) q[1];
rz(-pi) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.88586315) q[0];
sx q[0];
rz(-0.22561377) q[0];
sx q[0];
rz(-1.6560692) q[0];
rz(-pi) q[1];
rz(-0.26890484) q[2];
sx q[2];
rz(-0.78654002) q[2];
sx q[2];
rz(-1.5812909) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.4712217) q[1];
sx q[1];
rz(-2.2891217) q[1];
sx q[1];
rz(3.0918151) q[1];
x q[2];
rz(-0.66612996) q[3];
sx q[3];
rz(-2.1762848) q[3];
sx q[3];
rz(-0.44119409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3198118) q[2];
sx q[2];
rz(-2.767441) q[2];
sx q[2];
rz(-0.4846586) q[2];
rz(2.7759077) q[3];
sx q[3];
rz(-1.3644783) q[3];
sx q[3];
rz(1.1588089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0115688) q[0];
sx q[0];
rz(-0.31112177) q[0];
sx q[0];
rz(-0.097231641) q[0];
rz(0.10487996) q[1];
sx q[1];
rz(-2.361894) q[1];
sx q[1];
rz(1.3146776) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.275248) q[0];
sx q[0];
rz(-1.5731678) q[0];
sx q[0];
rz(-1.5771241) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.82523166) q[2];
sx q[2];
rz(-0.32352359) q[2];
sx q[2];
rz(-2.4458134) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-1.2995616) q[1];
sx q[1];
rz(-2.6202469) q[1];
sx q[1];
rz(-0.10592769) q[1];
rz(-pi) q[2];
rz(3.0544326) q[3];
sx q[3];
rz(-2.1059375) q[3];
sx q[3];
rz(-2.5453107) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(2.9647727) q[2];
sx q[2];
rz(-1.3975881) q[2];
sx q[2];
rz(-0.1532661) q[2];
rz(1.9249453) q[3];
sx q[3];
rz(-2.9235268) q[3];
sx q[3];
rz(-2.5254068) q[3];
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
rz(-2.0230873) q[0];
sx q[0];
rz(-3.1361134) q[0];
sx q[0];
rz(-1.5047005) q[0];
rz(-0.79832375) q[1];
sx q[1];
rz(-1.6183805) q[1];
sx q[1];
rz(-2.8062779) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17851795) q[0];
sx q[0];
rz(-1.8318684) q[0];
sx q[0];
rz(-0.065351323) q[0];
x q[1];
rz(2.3912001) q[2];
sx q[2];
rz(-1.5529035) q[2];
sx q[2];
rz(0.6070348) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.032051) q[1];
sx q[1];
rz(-0.45278835) q[1];
sx q[1];
rz(0.70346634) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.1481403) q[3];
sx q[3];
rz(-1.7029188) q[3];
sx q[3];
rz(0.41228774) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-3.0140784) q[2];
sx q[2];
rz(-1.1684343) q[2];
sx q[2];
rz(-0.43759313) q[2];
rz(-1.0420927) q[3];
sx q[3];
rz(-1.7362678) q[3];
sx q[3];
rz(2.5991345) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.9893148) q[0];
sx q[0];
rz(-0.88459009) q[0];
sx q[0];
rz(-2.6919795) q[0];
rz(-2.1647029) q[1];
sx q[1];
rz(-1.6376817) q[1];
sx q[1];
rz(0.50122112) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.86574449) q[0];
sx q[0];
rz(-1.8536708) q[0];
sx q[0];
rz(-0.19794835) q[0];
rz(-pi) q[1];
x q[1];
rz(0.21542549) q[2];
sx q[2];
rz(-1.1555724) q[2];
sx q[2];
rz(-1.2895467) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(0.93443643) q[1];
sx q[1];
rz(-1.109135) q[1];
sx q[1];
rz(-2.8769879) q[1];
x q[2];
rz(-1.1057749) q[3];
sx q[3];
rz(-1.2622202) q[3];
sx q[3];
rz(0.082482626) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-2.1180798) q[2];
sx q[2];
rz(-2.2492275) q[2];
sx q[2];
rz(-2.9238713) q[2];
rz(1.9067541) q[3];
sx q[3];
rz(-2.4775938) q[3];
sx q[3];
rz(-2.2209404) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49391838) q[0];
sx q[0];
rz(-1.4287345) q[0];
sx q[0];
rz(-2.7609974) q[0];
rz(-2.4299798) q[1];
sx q[1];
rz(-1.8743534) q[1];
sx q[1];
rz(1.4030917) q[1];
rz(3.0416476) q[2];
sx q[2];
rz(-0.47158416) q[2];
sx q[2];
rz(-0.33000962) q[2];
rz(-2.7357581) q[3];
sx q[3];
rz(-1.3530227) q[3];
sx q[3];
rz(1.4128662) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
