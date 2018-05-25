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
  <div style="width:40%;float:right"><img src="
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
