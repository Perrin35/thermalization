OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(-0.10652868) q[0];
sx q[0];
rz(-1.0892692) q[0];
sx q[0];
rz(0.16103345) q[0];
rz(1.610202) q[1];
sx q[1];
rz(-0.47709563) q[1];
sx q[1];
rz(0.49638003) q[1];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.30713233) q[0];
sx q[0];
rz(-1.4559329) q[0];
sx q[0];
rz(2.1644724) q[0];
rz(-pi) q[1];
x q[1];
rz(-1.7726937) q[2];
sx q[2];
rz(-1.9041833) q[2];
sx q[2];
rz(2.8436529) q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-1.5254933) q[1];
sx q[1];
rz(-1.4800737) q[1];
sx q[1];
rz(-2.4331122) q[1];
rz(-pi) q[2];
rz(-1.4708038) q[3];
sx q[3];
rz(-1.1649719) q[3];
sx q[3];
rz(2.0506746) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.51497841) q[2];
sx q[2];
rz(-0.86003059) q[2];
sx q[2];
rz(2.853945) q[2];
rz(1.3927762) q[3];
sx q[3];
rz(-0.98373047) q[3];
sx q[3];
rz(2.0387409) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.87067938) q[0];
sx q[0];
rz(-2.3864855) q[0];
sx q[0];
rz(0.57759181) q[0];
rz(-1.6060991) q[1];
sx q[1];
rz(-1.1178733) q[1];
sx q[1];
rz(1.5637406) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0751511) q[0];
sx q[0];
rz(-1.4581212) q[0];
sx q[0];
rz(0.22002797) q[0];
x q[1];
rz(-1.2969442) q[2];
sx q[2];
rz(-1.2784064) q[2];
sx q[2];
rz(-1.3016303) q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(2.8375081) q[1];
sx q[1];
rz(-1.4834852) q[1];
sx q[1];
rz(-1.1156032) q[1];
rz(-pi) q[2];
x q[2];
rz(-0.15474774) q[3];
sx q[3];
rz(-1.7259211) q[3];
sx q[3];
rz(2.250092) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.10721283) q[2];
sx q[2];
rz(-0.97637525) q[2];
sx q[2];
rz(1.0822901) q[2];
rz(-1.1632495) q[3];
sx q[3];
rz(-1.2000368) q[3];
sx q[3];
rz(0.64003402) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.0619693) q[0];
sx q[0];
rz(-2.847023) q[0];
sx q[0];
rz(2.9586155) q[0];
rz(3.0984763) q[1];
sx q[1];
rz(-0.95298195) q[1];
sx q[1];
rz(-1.9972237) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.17830081) q[0];
sx q[0];
rz(-2.460647) q[0];
sx q[0];
rz(2.2292024) q[0];
x q[1];
rz(-1.3470596) q[2];
sx q[2];
rz(-1.1615331) q[2];
sx q[2];
rz(0.70659107) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(pi/2) q[0];
rz(-2.548) q[1];
sx q[1];
rz(-1.9152194) q[1];
sx q[1];
rz(0.48732948) q[1];
rz(-1.5738437) q[3];
sx q[3];
rz(-1.2831732) q[3];
sx q[3];
rz(1.0685208) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-3.1220876) q[2];
sx q[2];
rz(-1.1818037) q[2];
sx q[2];
rz(-0.90467492) q[2];
rz(2.4217862) q[3];
sx q[3];
rz(-2.7147839) q[3];
sx q[3];
rz(-0.75511801) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
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
rz(1.4633789) q[0];
sx q[0];
rz(-2.1713874) q[0];
sx q[0];
rz(2.4257207) q[0];
rz(-0.32575682) q[1];
sx q[1];
rz(-1.0595067) q[1];
sx q[1];
rz(2.148927) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.0621322) q[0];
sx q[0];
rz(-1.8347164) q[0];
sx q[0];
rz(-0.53702766) q[0];
rz(-pi) q[1];
x q[1];
rz(-2.2797212) q[2];
sx q[2];
rz(-2.3700691) q[2];
sx q[2];
rz(-0.18415235) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-3.0992972) q[1];
sx q[1];
rz(-2.9534833) q[1];
sx q[1];
rz(2.6204965) q[1];
rz(2.643814) q[3];
sx q[3];
rz(-2.0911502) q[3];
sx q[3];
rz(-1.2342412) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(2.0585534) q[2];
sx q[2];
rz(-2.4251067) q[2];
sx q[2];
rz(-1.4286208) q[2];
rz(-2.2955017) q[3];
sx q[3];
rz(-1.6276136) q[3];
sx q[3];
rz(-1.2231474) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(pi/2) q[0];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.5457526) q[0];
sx q[0];
rz(-1.4056453) q[0];
sx q[0];
rz(2.3748421) q[0];
rz(-1.9013566) q[1];
sx q[1];
rz(-0.84754544) q[1];
sx q[1];
rz(0.3516745) q[1];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.2857916) q[0];
sx q[0];
rz(-1.8579146) q[0];
sx q[0];
rz(-0.65335269) q[0];
rz(-pi) q[1];
rz(-2.2064477) q[2];
sx q[2];
rz(-1.9872553) q[2];
sx q[2];
rz(0.93174975) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
x q[0];
rz(-2.5778351) q[1];
sx q[1];
rz(-1.1760684) q[1];
sx q[1];
rz(1.5549591) q[1];
rz(-2.8204927) q[3];
sx q[3];
rz(-1.0874815) q[3];
sx q[3];
rz(-1.4534284) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
rz(-2.9436283) q[2];
sx q[2];
rz(-1.4732271) q[2];
sx q[2];
rz(-2.6082805) q[2];
rz(-2.7126281) q[3];
sx q[3];
rz(-2.6140489) q[3];
sx q[3];
rz(-2.0146446) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi) q[1];
x q[2];
rz(-pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.11739843) q[0];
sx q[0];
rz(-2.0485931) q[0];
sx q[0];
rz(-0.54164106) q[0];
rz(-0.60846865) q[1];
sx q[1];
rz(-2.922373) q[1];
sx q[1];
rz(-2.0419962) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.8894419) q[0];
sx q[0];
rz(-0.60177207) q[0];
sx q[0];
rz(-0.88366951) q[0];
x q[1];
rz(2.9526688) q[2];
sx q[2];
rz(-1.0475698) q[2];
sx q[2];
rz(0.97666937) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(2.3179143) q[1];
sx q[1];
rz(-2.4412529) q[1];
sx q[1];
rz(2.8631696) q[1];
rz(2.3533456) q[3];
sx q[3];
rz(-1.8133834) q[3];
sx q[3];
rz(1.5188252) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.5114484) q[2];
sx q[2];
rz(-1.6161329) q[2];
sx q[2];
rz(2.8586094) q[2];
rz(2.4354368) q[3];
sx q[3];
rz(-1.599584) q[3];
sx q[3];
rz(-2.506536) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.054984897) q[0];
sx q[0];
rz(-0.98751706) q[0];
sx q[0];
rz(1.5299861) q[0];
rz(2.4967172) q[1];
sx q[1];
rz(-0.49352831) q[1];
sx q[1];
rz(2.9842916) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.29042127) q[0];
sx q[0];
rz(-0.84050814) q[0];
sx q[0];
rz(1.69676) q[0];
rz(-2.925161) q[2];
sx q[2];
rz(-2.2955403) q[2];
sx q[2];
rz(0.12491465) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(1.8017756) q[1];
sx q[1];
rz(-0.96313699) q[1];
sx q[1];
rz(-2.9031309) q[1];
x q[2];
rz(-0.39237202) q[3];
sx q[3];
rz(-2.4761768) q[3];
sx q[3];
rz(0.80442807) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-1.1823696) q[2];
sx q[2];
rz(-2.5490641) q[2];
sx q[2];
rz(0.83089337) q[2];
rz(-0.043878555) q[3];
sx q[3];
rz(-1.2336122) q[3];
sx q[3];
rz(-0.20251814) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(pi/2) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.6458994) q[0];
sx q[0];
rz(-0.93911397) q[0];
sx q[0];
rz(3.1307401) q[0];
rz(2.8649578) q[1];
sx q[1];
rz(-2.3697772) q[1];
sx q[1];
rz(2.8483134) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5787443) q[0];
sx q[0];
rz(-2.1132073) q[0];
sx q[0];
rz(-1.5777274) q[0];
rz(-pi) q[1];
rz(2.2145055) q[2];
sx q[2];
rz(-1.7037399) q[2];
sx q[2];
rz(-1.5333652) q[2];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi/2) q[0];
rz(-1.7640904) q[1];
sx q[1];
rz(-1.1524156) q[1];
sx q[1];
rz(2.3368895) q[1];
rz(0.46383143) q[3];
sx q[3];
rz(-1.9411191) q[3];
sx q[3];
rz(0.90255373) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(3.098112) q[2];
sx q[2];
rz(-0.094823368) q[2];
sx q[2];
rz(-3.0528255) q[2];
rz(-0.25012112) q[3];
sx q[3];
rz(-1.1238031) q[3];
sx q[3];
rz(0.5789825) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
sx q[3];
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
rz(-1.0042689) q[0];
sx q[0];
rz(-2.9172638) q[0];
sx q[0];
rz(-2.0776757) q[0];
rz(2.6990199) q[1];
sx q[1];
rz(-1.4241709) q[1];
sx q[1];
rz(-1.7907422) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.0161184) q[0];
sx q[0];
rz(-0.11575143) q[0];
sx q[0];
rz(-2.3257757) q[0];
x q[1];
rz(2.4814018) q[2];
sx q[2];
rz(-1.5056416) q[2];
sx q[2];
rz(-2.8934663) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(-2.0296214) q[1];
sx q[1];
rz(-2.0051149) q[1];
sx q[1];
rz(3.1236468) q[1];
rz(-pi) q[2];
rz(2.001858) q[3];
sx q[3];
rz(-0.24340478) q[3];
sx q[3];
rz(0.36153015) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-2.0315447) q[2];
sx q[2];
rz(-1.8507379) q[2];
sx q[2];
rz(-0.44719493) q[2];
rz(1.7556919) q[3];
sx q[3];
rz(-2.8409676) q[3];
sx q[3];
rz(0.51469222) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(pi/2) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.31496012) q[0];
sx q[0];
rz(-1.6199912) q[0];
sx q[0];
rz(0.78053027) q[0];
rz(-0.42778095) q[1];
sx q[1];
rz(-0.3573187) q[1];
sx q[1];
rz(0.16960493) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.0376301) q[0];
sx q[0];
rz(-2.5314405) q[0];
sx q[0];
rz(2.2274341) q[0];
rz(-0.054762997) q[2];
sx q[2];
rz(-2.7191396) q[2];
sx q[2];
rz(-1.7516608) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.5055713) q[1];
sx q[1];
rz(-0.95658703) q[1];
sx q[1];
rz(-3.0843656) q[1];
rz(-pi) q[2];
rz(-1.2492368) q[3];
sx q[3];
rz(-1.4351298) q[3];
sx q[3];
rz(-2.3615169) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.2910989) q[2];
sx q[2];
rz(-1.9591745) q[2];
sx q[2];
rz(0.69236857) q[2];
rz(0.46323562) q[3];
sx q[3];
rz(-1.654947) q[3];
sx q[3];
rz(1.108981) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
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
rz(pi/2) q[0];
x q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.8511843) q[0];
sx q[0];
rz(-1.52607) q[0];
sx q[0];
rz(-2.4551256) q[0];
rz(0.51207536) q[1];
sx q[1];
rz(-0.33437406) q[1];
sx q[1];
rz(-2.1127111) q[1];
rz(2.8725273) q[2];
sx q[2];
rz(-1.27956) q[2];
sx q[2];
rz(0.051824311) q[2];
rz(-1.7789755) q[3];
sx q[3];
rz(-1.4559742) q[3];
sx q[3];
rz(-2.5696587) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];
