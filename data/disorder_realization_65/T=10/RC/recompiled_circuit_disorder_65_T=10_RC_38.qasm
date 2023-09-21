OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-2.8310175) q[0];
sx q[0];
rz(-2.0895045) q[0];
sx q[0];
rz(-1.6488099) q[0];
rz(-2.2812023) q[1];
sx q[1];
rz(-0.69505039) q[1];
sx q[1];
rz(1.0746497) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2225392) q[0];
sx q[0];
rz(-3.1226282) q[0];
sx q[0];
rz(2.8595964) q[0];
rz(-pi) q[1];
x q[1];
rz(1.7311814) q[2];
sx q[2];
rz(-0.78750247) q[2];
sx q[2];
rz(0.10057848) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.5584471) q[1];
sx q[1];
rz(-1.4902643) q[1];
sx q[1];
rz(-3.0196378) q[1];
rz(2.4514276) q[3];
sx q[3];
rz(-0.65342045) q[3];
sx q[3];
rz(2.0131468) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(2.3124714) q[2];
sx q[2];
rz(-2.1437936) q[2];
sx q[2];
rz(-1.6281698) q[2];
rz(-1.9681905) q[3];
sx q[3];
rz(-1.6956001) q[3];
sx q[3];
rz(-2.9287958) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5728977) q[0];
sx q[0];
rz(-2.499489) q[0];
sx q[0];
rz(-1.2233541) q[0];
rz(-2.9808295) q[1];
sx q[1];
rz(-1.5784966) q[1];
sx q[1];
rz(1.6404023) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.879724) q[0];
sx q[0];
rz(-0.12517087) q[0];
sx q[0];
rz(-0.20010389) q[0];
rz(2.03962) q[2];
sx q[2];
rz(-1.2906244) q[2];
sx q[2];
rz(0.75129347) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(1.3850296) q[1];
sx q[1];
rz(-0.90365138) q[1];
sx q[1];
rz(-1.7900311) q[1];
rz(-pi) q[2];
x q[2];
rz(-1.1797446) q[3];
sx q[3];
rz(-1.565233) q[3];
sx q[3];
rz(0.55913505) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(-1.8555277) q[2];
sx q[2];
rz(-2.3748886) q[2];
sx q[2];
rz(-1.6274874) q[2];
rz(1.9747915) q[3];
sx q[3];
rz(-1.5224001) q[3];
sx q[3];
rz(-2.2272026) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0062155) q[0];
sx q[0];
rz(-1.051798) q[0];
sx q[0];
rz(1.0789543) q[0];
rz(-1.1795993) q[1];
sx q[1];
rz(-2.1005354) q[1];
sx q[1];
rz(-1.5037781) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30341879) q[0];
sx q[0];
rz(-1.5056599) q[0];
sx q[0];
rz(1.6584048) q[0];
rz(-pi) q[1];
x q[1];
rz(2.1554699) q[2];
sx q[2];
rz(-1.976527) q[2];
sx q[2];
rz(-0.090678064) q[2];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(1.9864195) q[1];
sx q[1];
rz(-1.5281614) q[1];
sx q[1];
rz(-2.8405872) q[1];
x q[2];
rz(0.13111968) q[3];
sx q[3];
rz(-0.91851202) q[3];
sx q[3];
rz(2.2281447) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(0.43061313) q[2];
sx q[2];
rz(-2.6617699) q[2];
sx q[2];
rz(-0.7437931) q[2];
rz(2.8841694) q[3];
sx q[3];
rz(-2.0580723) q[3];
sx q[3];
rz(-0.58810365) q[3];
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
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.96823111) q[0];
sx q[0];
rz(-3.0259202) q[0];
sx q[0];
rz(-1.051735) q[0];
rz(2.5412718) q[1];
sx q[1];
rz(-0.61906639) q[1];
sx q[1];
rz(-1.1490885) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.34654348) q[0];
sx q[0];
rz(-1.330089) q[0];
sx q[0];
rz(1.4562777) q[0];
rz(-pi) q[1];
rz(-0.97781424) q[2];
sx q[2];
rz(-1.2056418) q[2];
sx q[2];
rz(-2.1821373) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-0.035605343) q[1];
sx q[1];
rz(-2.3590901) q[1];
sx q[1];
rz(-1.993202) q[1];
rz(-pi) q[2];
rz(1.2265008) q[3];
sx q[3];
rz(-2.2831884) q[3];
sx q[3];
rz(-0.32574367) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-1.8950618) q[2];
sx q[2];
rz(-1.6677413) q[2];
sx q[2];
rz(2.4728298) q[2];
rz(-1.8956005) q[3];
sx q[3];
rz(-0.99530882) q[3];
sx q[3];
rz(-3.1252089) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
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
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.65288654) q[0];
sx q[0];
rz(-1.1504494) q[0];
sx q[0];
rz(-2.0892129) q[0];
rz(1.2106238) q[1];
sx q[1];
rz(-1.4167507) q[1];
sx q[1];
rz(-2.2573684) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0563755) q[0];
sx q[0];
rz(-1.7612805) q[0];
sx q[0];
rz(-0.18116829) q[0];
x q[1];
rz(0.82489478) q[2];
sx q[2];
rz(-2.5755304) q[2];
sx q[2];
rz(2.538946) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(2.4864038) q[1];
sx q[1];
rz(-0.15895325) q[1];
sx q[1];
rz(2.7093191) q[1];
x q[2];
rz(-2.103881) q[3];
sx q[3];
rz(-1.248705) q[3];
sx q[3];
rz(1.0853634) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(0.23652442) q[2];
sx q[2];
rz(-0.69513598) q[2];
sx q[2];
rz(2.0489676) q[2];
rz(0.24400273) q[3];
sx q[3];
rz(-2.2677393) q[3];
sx q[3];
rz(-0.75105914) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[0];
x q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.4051751) q[0];
sx q[0];
rz(-0.002679499) q[0];
sx q[0];
rz(1.3177692) q[0];
rz(-2.4353943) q[1];
sx q[1];
rz(-1.6786007) q[1];
sx q[1];
rz(2.172519) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.58020619) q[0];
sx q[0];
rz(-1.5302587) q[0];
sx q[0];
rz(1.4617306) q[0];
rz(-pi) q[1];
rz(2.4275052) q[2];
sx q[2];
rz(-1.6220777) q[2];
sx q[2];
rz(-0.72061348) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.6094738) q[1];
sx q[1];
rz(-1.1540968) q[1];
sx q[1];
rz(-2.6344224) q[1];
rz(-3*pi/11) q[3];
sx q[3];
rz(-0.65398765) q[3];
sx q[3];
rz(1.4531144) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.6270854) q[2];
sx q[2];
rz(-0.58834499) q[2];
sx q[2];
rz(-2.7039841) q[2];
rz(0.058241025) q[3];
sx q[3];
rz(-2.2831459) q[3];
sx q[3];
rz(1.6102012) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3095734) q[0];
sx q[0];
rz(-1.0252527) q[0];
sx q[0];
rz(1.3597885) q[0];
rz(-1.6925905) q[1];
sx q[1];
rz(-0.89869181) q[1];
sx q[1];
rz(1.4356027) q[1];
rz(-pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.3754417) q[0];
sx q[0];
rz(-1.0784082) q[0];
sx q[0];
rz(2.173461) q[0];
rz(-pi) q[1];
x q[1];
rz(0.82201634) q[2];
sx q[2];
rz(-0.026958131) q[2];
sx q[2];
rz(-0.4617304) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(0.85305271) q[1];
sx q[1];
rz(-1.2097881) q[1];
sx q[1];
rz(0.32220528) q[1];
x q[2];
rz(-2.0766047) q[3];
sx q[3];
rz(-1.1899663) q[3];
sx q[3];
rz(1.8246458) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[1];
rz(-0.74063611) q[2];
sx q[2];
rz(-1.8813958) q[2];
sx q[2];
rz(2.9675972) q[2];
rz(1.0081572) q[3];
sx q[3];
rz(-0.62767902) q[3];
sx q[3];
rz(-0.15667668) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.85132861) q[0];
sx q[0];
rz(-0.86306089) q[0];
sx q[0];
rz(3.0086349) q[0];
rz(0.48775396) q[1];
sx q[1];
rz(-2.1552591) q[1];
sx q[1];
rz(1.3652323) q[1];
rz(pi/2) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4505969) q[0];
sx q[0];
rz(-2.9857675) q[0];
sx q[0];
rz(1.063892) q[0];
rz(1.7403142) q[2];
sx q[2];
rz(-3.0119544) q[2];
sx q[2];
rz(0.6946176) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-1.8147162) q[1];
sx q[1];
rz(-1.6232345) q[1];
sx q[1];
rz(0.53175064) q[1];
rz(-pi) q[2];
x q[2];
rz(-2.7435999) q[3];
sx q[3];
rz(-2.6905305) q[3];
sx q[3];
rz(-3.1068902) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(0.075295538) q[2];
sx q[2];
rz(-1.8981045) q[2];
sx q[2];
rz(0.42262849) q[2];
rz(2.1832809) q[3];
sx q[3];
rz(-0.13923968) q[3];
sx q[3];
rz(-2.9039834) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.002710297) q[0];
sx q[0];
rz(-2.8399828) q[0];
sx q[0];
rz(2.8310006) q[0];
rz(0.25767913) q[1];
sx q[1];
rz(-1.5035166) q[1];
sx q[1];
rz(0.92393595) q[1];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.63147488) q[0];
sx q[0];
rz(-2.3259813) q[0];
sx q[0];
rz(-0.60879137) q[0];
rz(-pi) q[1];
rz(3.1191101) q[2];
sx q[2];
rz(-0.48869952) q[2];
sx q[2];
rz(1.894001) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.9388705) q[1];
sx q[1];
rz(-2.466423) q[1];
sx q[1];
rz(-0.90941456) q[1];
rz(-pi) q[2];
rz(-0.71420788) q[3];
sx q[3];
rz(-2.0651157) q[3];
sx q[3];
rz(-0.89356542) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.6442287) q[2];
sx q[2];
rz(-0.53933898) q[2];
sx q[2];
rz(1.3941992) q[2];
rz(1.9723069) q[3];
sx q[3];
rz(-1.0401657) q[3];
sx q[3];
rz(-2.6031301) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
sx q[3];
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
rz(0.37344638) q[0];
sx q[0];
rz(-3.0420619) q[0];
sx q[0];
rz(-1.753153) q[0];
rz(0.0069847981) q[1];
sx q[1];
rz(-0.81022898) q[1];
sx q[1];
rz(3.1034234) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.5157493) q[0];
sx q[0];
rz(-2.1793097) q[0];
sx q[0];
rz(0.032511062) q[0];
rz(-pi) q[1];
x q[1];
rz(-0.091911749) q[2];
sx q[2];
rz(-1.1895864) q[2];
sx q[2];
rz(3.0222169) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(0.2529777) q[1];
sx q[1];
rz(-1.9441609) q[1];
sx q[1];
rz(-0.044815973) q[1];
x q[2];
rz(2.2977703) q[3];
sx q[3];
rz(-2.3355964) q[3];
sx q[3];
rz(-0.7741011) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(0.96597153) q[2];
sx q[2];
rz(-1.1493382) q[2];
sx q[2];
rz(1.9434628) q[2];
rz(1.0839869) q[3];
sx q[3];
rz(-0.67892781) q[3];
sx q[3];
rz(-2.9057124) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[2];
x q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.1643628) q[0];
sx q[0];
rz(-1.2156163) q[0];
sx q[0];
rz(-1.6396133) q[0];
rz(-1.4981131) q[1];
sx q[1];
rz(-2.5479981) q[1];
sx q[1];
rz(2.610276) q[1];
rz(0.48508587) q[2];
sx q[2];
rz(-3.017364) q[2];
sx q[2];
rz(-0.45807522) q[2];
rz(2.8244143) q[3];
sx q[3];
rz(-0.89430292) q[3];
sx q[3];
rz(3.111307) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];