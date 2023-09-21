OPENQASM 2.0;
include "qelib1.inc";
qreg q[4];
creg c[4];
sx q[0];
sx q[1];
sx q[2];
sx q[3];
barrier q[0],q[1],q[2],q[3];
rz(2.3611276) q[0];
sx q[0];
rz(-0.68929231) q[0];
sx q[0];
rz(2.8110992) q[0];
rz(-2.8117872) q[1];
sx q[1];
rz(-2.2916315) q[1];
sx q[1];
rz(2.4324774) q[1];
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.56624428) q[0];
sx q[0];
rz(-0.69501221) q[0];
sx q[0];
rz(-2.038875) q[0];
rz(-pi) q[1];
x q[1];
rz(1.2572631) q[2];
sx q[2];
rz(-1.5402113) q[2];
sx q[2];
rz(2.7262296) q[2];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(-1.658537) q[1];
sx q[1];
rz(-0.72209789) q[1];
sx q[1];
rz(-2.2621821) q[1];
rz(-pi) q[2];
rz(0.37681864) q[3];
sx q[3];
rz(-1.0816649) q[3];
sx q[3];
rz(2.7917142) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi/2) q[1];
rz(1.7314529) q[2];
sx q[2];
rz(-2.7316766) q[2];
sx q[2];
rz(-1.5343792) q[2];
rz(0.93506995) q[3];
sx q[3];
rz(-1.3053223) q[3];
sx q[3];
rz(0.7888166) q[3];
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
x q[0];
rz(pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.4193029) q[0];
sx q[0];
rz(-3.0794444) q[0];
sx q[0];
rz(0.62227917) q[0];
rz(0.17624804) q[1];
sx q[1];
rz(-1.2146249) q[1];
sx q[1];
rz(-0.91631779) q[1];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.5641862) q[0];
sx q[0];
rz(-0.86125492) q[0];
sx q[0];
rz(-0.90155154) q[0];
x q[1];
rz(1.4257405) q[2];
sx q[2];
rz(-2.3999891) q[2];
sx q[2];
rz(1.8650101) q[2];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(-1.2434477) q[1];
sx q[1];
rz(-1.2444278) q[1];
sx q[1];
rz(1.9279724) q[1];
x q[2];
rz(1.6244435) q[3];
sx q[3];
rz(-2.7430153) q[3];
sx q[3];
rz(0.22565354) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-2.9006485) q[2];
sx q[2];
rz(-0.99664656) q[2];
sx q[2];
rz(-0.82143482) q[2];
rz(-0.017283043) q[3];
sx q[3];
rz(-0.79764962) q[3];
sx q[3];
rz(0.78330529) q[3];
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
rz(pi/2) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-pi) q[2];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.18773742) q[0];
sx q[0];
rz(-1.563235) q[0];
sx q[0];
rz(1.0082555) q[0];
rz(-3.1058274) q[1];
sx q[1];
rz(-1.6556031) q[1];
sx q[1];
rz(0.52454138) q[1];
rz(-pi/2) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.6150104) q[0];
sx q[0];
rz(-2.9752762) q[0];
sx q[0];
rz(1.7517356) q[0];
x q[1];
rz(-0.67851615) q[2];
sx q[2];
rz(-1.3635474) q[2];
sx q[2];
rz(-3.0845272) q[2];
rz(-pi/2) q[3];
sx q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-0.14568612) q[1];
sx q[1];
rz(-1.5090764) q[1];
sx q[1];
rz(1.5476336) q[1];
x q[2];
rz(-1.8936162) q[3];
sx q[3];
rz(-0.4399235) q[3];
sx q[3];
rz(-2.3390714) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(2.3699469) q[2];
sx q[2];
rz(-1.5180072) q[2];
sx q[2];
rz(1.6463722) q[2];
rz(1.8203991) q[3];
sx q[3];
rz(-1.7318055) q[3];
sx q[3];
rz(-0.43740073) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(-pi) q[1];
rz(-pi) q[2];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(-pi) q[1];
x q[1];
rz(pi/2) q[2];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8111073) q[0];
sx q[0];
rz(-2.2225668) q[0];
sx q[0];
rz(-3.0526429) q[0];
rz(-0.51070172) q[1];
sx q[1];
rz(-2.316244) q[1];
sx q[1];
rz(-0.68960062) q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3084429) q[0];
sx q[0];
rz(-1.8652548) q[0];
sx q[0];
rz(1.4008646) q[0];
x q[1];
rz(0.34822779) q[2];
sx q[2];
rz(-0.59154445) q[2];
sx q[2];
rz(2.0915742) q[2];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi/2) q[0];
rz(-1.2526306) q[1];
sx q[1];
rz(-1.8969715) q[1];
sx q[1];
rz(0.77146448) q[1];
rz(-pi) q[2];
rz(2.6570286) q[3];
sx q[3];
rz(-1.441941) q[3];
sx q[3];
rz(-1.3467005) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
x q[1];
rz(2.4758063) q[2];
sx q[2];
rz(-0.67725724) q[2];
sx q[2];
rz(-2.518667) q[2];
rz(-2.0056491) q[3];
sx q[3];
rz(-0.1853075) q[3];
sx q[3];
rz(-0.46666551) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
rz(-pi) q[1];
x q[1];
rz(-pi) q[2];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.85161197) q[0];
sx q[0];
rz(-1.2427793) q[0];
sx q[0];
rz(-0.583453) q[0];
rz(-1.9955697) q[1];
sx q[1];
rz(-1.4910411) q[1];
sx q[1];
rz(-1.4978283) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(1.0935055) q[0];
sx q[0];
rz(-1.830528) q[0];
sx q[0];
rz(-0.78608677) q[0];
rz(-pi) q[1];
x q[1];
rz(-3.0450902) q[2];
sx q[2];
rz(-0.8000024) q[2];
sx q[2];
rz(-0.9347136) q[2];
x q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
sx q[0];
rz(pi/2) q[0];
rz(1.076482) q[1];
sx q[1];
rz(-1.5981734) q[1];
sx q[1];
rz(-0.70365023) q[1];
rz(-pi) q[2];
rz(1.7329526) q[3];
sx q[3];
rz(-2.080653) q[3];
sx q[3];
rz(1.4354524) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(-0.65537611) q[2];
sx q[2];
rz(-1.5465753) q[2];
sx q[2];
rz(-1.5779457) q[2];
rz(-0.90562138) q[3];
sx q[3];
rz(-2.8712397) q[3];
sx q[3];
rz(-2.5031228) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
sx q[0];
rz(-pi) q[0];
sx q[2];
rz(-pi) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.49801302) q[0];
sx q[0];
rz(-1.6828515) q[0];
sx q[0];
rz(0.046982732) q[0];
rz(2.9934096) q[1];
sx q[1];
rz(-0.89887416) q[1];
sx q[1];
rz(1.7061589) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3008266) q[0];
sx q[0];
rz(-0.78221417) q[0];
sx q[0];
rz(-1.7813111) q[0];
rz(-0.058850364) q[2];
sx q[2];
rz(-2.119679) q[2];
sx q[2];
rz(2.5318052) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
x q[0];
rz(2.6378577) q[1];
sx q[1];
rz(-1.6956455) q[1];
sx q[1];
rz(-1.9810852) q[1];
rz(-pi) q[2];
x q[2];
rz(1.3330323) q[3];
sx q[3];
rz(-1.5441582) q[3];
sx q[3];
rz(0.8386855) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(pi/2) q[1];
rz(1.0753714) q[2];
sx q[2];
rz(-2.1263289) q[2];
sx q[2];
rz(0.071468778) q[2];
rz(1.6890769) q[3];
sx q[3];
rz(-1.6966597) q[3];
sx q[3];
rz(2.215109) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
x q[1];
rz(-pi/2) q[1];
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
x q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-2.8554095) q[0];
sx q[0];
rz(-1.8137285) q[0];
sx q[0];
rz(-2.563971) q[0];
rz(-1.2795992) q[1];
sx q[1];
rz(-1.3459233) q[1];
sx q[1];
rz(1.0095899) q[1];
rz(-pi) q[2];
x q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.6972835) q[0];
sx q[0];
rz(-2.5464006) q[0];
sx q[0];
rz(-0.62023456) q[0];
rz(-pi) q[1];
rz(-1.3426443) q[2];
sx q[2];
rz(-1.1953127) q[2];
sx q[2];
rz(2.0827039) q[2];
rz(pi/2) q[3];
sx q[3];
rz(-pi) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi) q[0];
rz(0.84570388) q[1];
sx q[1];
rz(-2.0659975) q[1];
sx q[1];
rz(1.1843029) q[1];
rz(-pi) q[2];
rz(0.61543492) q[3];
sx q[3];
rz(-2.2781567) q[3];
sx q[3];
rz(2.5497041) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(pi/2) q[1];
sx q[1];
rz(-2.0916831) q[2];
sx q[2];
rz(-2.78806) q[2];
sx q[2];
rz(-1.9419149) q[2];
rz(-2.6692634) q[3];
sx q[3];
rz(-1.6475369) q[3];
sx q[3];
rz(2.1388617) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(pi/2) q[1];
rz(-pi) q[2];
x q[2];
x q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
x q[0];
x q[1];
sx q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.49884477) q[0];
sx q[0];
rz(-1.6413178) q[0];
sx q[0];
rz(-0.2391267) q[0];
rz(-2.3587976) q[1];
sx q[1];
rz(-0.51819623) q[1];
sx q[1];
rz(2.696864) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(2.3702104) q[0];
sx q[0];
rz(-1.0078197) q[0];
sx q[0];
rz(2.6315106) q[0];
rz(-pi) q[1];
x q[1];
rz(0.93034805) q[2];
sx q[2];
rz(-2.1313416) q[2];
sx q[2];
rz(-2.408574) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-pi) q[0];
rz(-1.7263033) q[1];
sx q[1];
rz(-2.2518034) q[1];
sx q[1];
rz(0.16420941) q[1];
rz(-0.29856155) q[3];
sx q[3];
rz(-0.42793722) q[3];
sx q[3];
rz(3.0446133) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(-pi) q[1];
rz(0.42177054) q[2];
sx q[2];
rz(-0.82227102) q[2];
sx q[2];
rz(2.3262809) q[2];
rz(1.0347962) q[3];
sx q[3];
rz(-1.5650322) q[3];
sx q[3];
rz(-1.3172654) q[3];
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
x q[0];
rz(-pi/2) q[0];
x q[1];
x q[2];
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-1.3141044) q[0];
sx q[0];
rz(-1.0359456) q[0];
sx q[0];
rz(-2.8544193) q[0];
rz(-0.18889591) q[1];
sx q[1];
rz(-2.4177528) q[1];
sx q[1];
rz(-2.8093991) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.19776343) q[0];
sx q[0];
rz(-1.257886) q[0];
sx q[0];
rz(0.80425941) q[0];
rz(1.6261149) q[2];
sx q[2];
rz(-0.69960591) q[2];
sx q[2];
rz(2.827022) q[2];
rz(-pi) q[3];
sx q[3];
rz(pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(-pi/2) q[0];
sx q[0];
rz(-3.0041973) q[1];
sx q[1];
rz(-1.157837) q[1];
sx q[1];
rz(-1.2910299) q[1];
rz(1.9751777) q[3];
sx q[3];
rz(-1.0746733) q[3];
sx q[3];
rz(0.29508428) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi) q[1];
sx q[1];
rz(pi/2) q[1];
rz(2.2087848) q[2];
sx q[2];
rz(-2.1366426) q[2];
sx q[2];
rz(2.3020111) q[2];
rz(1.8509289) q[3];
sx q[3];
rz(-1.3495812) q[3];
sx q[3];
rz(0.65657842) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
x q[1];
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
rz(-pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.055450913) q[0];
sx q[0];
rz(-0.76776004) q[0];
sx q[0];
rz(-2.0822051) q[0];
rz(2.1620031) q[1];
sx q[1];
rz(-1.1971808) q[1];
sx q[1];
rz(-2.8870781) q[1];
sx q[2];
rz(pi/2) q[2];
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(0.33728889) q[0];
sx q[0];
rz(-2.8992607) q[0];
sx q[0];
rz(-1.8810349) q[0];
rz(-pi) q[1];
x q[1];
rz(2.6539831) q[2];
sx q[2];
rz(-1.7874103) q[2];
sx q[2];
rz(0.85607869) q[2];
rz(-pi) q[3];
sx q[3];
rz(-pi/2) q[3];
barrier q[0],q[3],q[2],q[1];
cx q[2],q[1];
barrier q[0],q[3],q[2],q[1];
rz(pi/2) q[0];
sx q[0];
rz(2.3855615) q[1];
sx q[1];
rz(-0.33547151) q[1];
sx q[1];
rz(2.8940593) q[1];
rz(-pi) q[2];
x q[2];
rz(2.7032095) q[3];
sx q[3];
rz(-1.66072) q[3];
sx q[3];
rz(0.59059483) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
rz(-pi/2) q[1];
sx q[1];
rz(2.09482) q[2];
sx q[2];
rz(-0.54868713) q[2];
sx q[2];
rz(2.0937031) q[2];
rz(1.3607599) q[3];
sx q[3];
rz(-2.4436617) q[3];
sx q[3];
rz(2.5361983) q[3];
barrier q[1],q[3],q[2];
cx q[3],q[2];
barrier q[1],q[3],q[2];
sx q[1];
rz(-pi/2) q[3];
sx q[3];
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
barrier q[2],q[1],q[0];
cx q[1],q[0];
barrier q[2],q[1],q[0];
rz(-0.027325252) q[0];
sx q[0];
rz(-1.1189168) q[0];
sx q[0];
rz(3.0615321) q[0];
rz(0.36021532) q[1];
sx q[1];
rz(-1.6811014) q[1];
sx q[1];
rz(-0.95492687) q[1];
rz(0.24047492) q[2];
sx q[2];
rz(-2.0839811) q[2];
sx q[2];
rz(2.6032084) q[2];
rz(2.9110254) q[3];
sx q[3];
rz(-2.1721526) q[3];
sx q[3];
rz(-0.9823907) q[3];
barrier q[0],q[1],q[2],q[3];
measure q[0] -> c[0];
measure q[1] -> c[1];
measure q[2] -> c[2];
measure q[3] -> c[3];