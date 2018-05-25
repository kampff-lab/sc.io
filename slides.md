---
layout: slides
title: Slides
permalink: /slides/
---

<section>
<h1>SC.IO</h1>
<h2>Scientific Input-Output;<br/>Science-Open</h2>
<a href="http://www.ucl.ac.uk/swc">
  <img style="width:50%" class="plain" src="{{ site.baseurl }}/assets/images/swc.png">
</a>
</section>

{% for slide in site.intro %}
<section>
<h3>{{ slide.title }}</h3>
{% if slide.image %}
<div style="width:100%">
  <div style="width:60%;float:left">{{ slide.content }}</div>
  <div style="width:40%;float:right"><img class="plain" src="
{% if slide.image contains '://' %}
{{slide.image}}
{% else %}
{{ slide.image | prepend: site.baseurl }}
{% endif %}" />
  </div>
</div>
{% else %}
{{ slide.content }}
{% endif %}
</section>
{% endfor %}

<section>
  <section>
  <h1>Round-tables</h1>
  <h5 class="fragment">Requirements</h5>
  <h5 class="fragment">Experiences</h5>
  <h5 class="fragment">Obstacles</h5>
  </section>

  <section>
  <h2>Manuscripts as living documents</h2>
  <h5></h5>
  </section>

  <section>
  <h2>Continuous, transparent pre- and peer-review</h2>
  <h5></h5>
  </section>

  <section>
  <h2>Hosting and validating raw data</h2>
  <h5></h5>
  </section>

  <section>
  <h2>Reproducibility of analysis code and intermediate outputs</h2>
  <h5></h5>
  </section>

  <section>
  <h2>Strategies for reporting non-digital contributions</h2>
  <h5></h5>
  </section>

  <section>
  <h2>Metrics for evaluating impact of open-science methodologies</h2>
  <h5></h5>
  </section>
</section>

<section>
  <section>
  <h1>Open-ended discussions</h1>
  </section>

  <section>
  <h2>Scientific ecosystem in an open-science world</h2>
  </section>

  <section>
  <h2>Negative and exploratory results</h2>
  </section>

  <section>
  <h2>Applying open-science philosophy to scientific social infrastructure</h2>
  </section>

  <section>
  <h2>Increasing accessibility of open-science tools</h2>
  </section>

  <section>
  <h2>Wrap-up / Next steps</h2>
  </section>
</section>